"""Safe argv-based external process execution with immutable JSON records."""
from __future__ import annotations
from dataclasses import asdict, dataclass
from enum import Enum
import os, shlex, signal, subprocess, threading, time
from pathlib import Path
from typing import Mapping
from .errors import CommandExecutionError, CommandTimeoutError, ExpectedOutputMissingError
from .provenance import atomic_write_json, file_identity, utc_now
from .tools import ContainerMetadata, ExecutionMode, Tool
from .version import pipeline_identity

class ExecutionStatus(str, Enum):
    PLANNED="PLANNED"; RUNNING="RUNNING"; SUCCEEDED="SUCCEEDED"; FAILED="FAILED"; TIMED_OUT="TIMED_OUT"; CANCELLED="CANCELLED"; DRY_RUN="DRY_RUN"; OUTPUT_VALIDATION_FAILED="OUTPUT_VALIDATION_FAILED"
@dataclass(frozen=True)
class InputDeclaration: input_id:str; path:str; required:bool=True; role:str|None=None
@dataclass(frozen=True)
class OutputDeclaration: output_id:str; path:str; required:bool=True; allow_empty:bool=False; media_type:str|None=None; role:str|None=None
@dataclass(frozen=True)
class CommandSpec:
    command_id:str; stage_id:str; tool_id:str; argv:tuple[str,...]; working_directory:str
    environment_overrides:Mapping[str,str]=None; redacted_environment_keys:tuple[str,...]=(); declared_inputs:tuple[InputDeclaration,...]=(); declared_outputs:tuple[OutputDeclaration,...]=()
    expected_exit_codes:tuple[int,...]=(0,); timeout_seconds:float|None=None; overwrite:bool=False; execution_mode:ExecutionMode=ExecutionMode.NATIVE; container_metadata:ContainerMetadata|None=None
    def __post_init__(self):
        object.__setattr__(self,"environment_overrides",dict(self.environment_overrides or {}))
        if not Path(self.working_directory).is_absolute(): raise ValueError("working_directory must be absolute")

class CancellationToken:
    def __init__(self): self._event=threading.Event()
    def cancel(self): self._event.set()
    @property
    def cancelled(self): return self._event.is_set()

def _serialize(items): return [asdict(x) for x in items]
def execute(spec:CommandSpec, tool:Tool, record_path:str|Path, log_root:str|Path, *, dry_run=False, cancellation:CancellationToken|None=None, termination_grace_seconds=.5):
    work=Path(spec.working_directory); work.mkdir(parents=True,exist_ok=True)
    log_dir=Path(log_root)/spec.stage_id; log_dir.mkdir(parents=True,exist_ok=True)
    stdout_path=log_dir/f"{spec.command_id}.stdout.log"; stderr_path=log_dir/f"{spec.command_id}.stderr.log"; combined_path=log_dir/f"{spec.command_id}.combined.log"
    argv=spec.argv
    container_identity=None
    if spec.execution_mode is ExecutionMode.APPTAINER:
        if not spec.container_metadata: raise ValueError("APPTAINER requires container_metadata")
        argv=spec.container_metadata.host_argv(spec.argv,spec.working_directory)
        container_identity={"container_image":spec.container_metadata.container_image,"container_sha256":spec.container_metadata.container_sha256,"bind_mounts":list(spec.container_metadata.bind_mounts),"internal_executable":spec.container_metadata.internal_executable,"internal_argv":list(spec.argv),"host_argv":list(argv),"host_working_directory":spec.working_directory,"container_working_directory":spec.container_metadata.container_working_directory}
    existing=[o.path for o in spec.declared_outputs if Path(o.path).exists()]
    if existing and not spec.overwrite: raise CommandExecutionError("declared output already exists",existing_outputs=existing)
    identity=pipeline_identity(); started=utc_now(); begin=time.monotonic()
    env_visible={k:("<REDACTED>" if k in spec.redacted_environment_keys else v) for k,v in spec.environment_overrides.items()}
    base={"record_schema_version":"1.0","command_id":spec.command_id,"stage_id":spec.stage_id,"tool_id":spec.tool_id,"tool_display_name":tool.display_name,"configured_executable":tool.configured_executable,"resolved_executable":tool.resolved_executable,"tool_version":tool.detected_version,"pipeline_version":identity.pipeline_version,"pipeline_git_commit":identity.git_commit,"pipeline_git_dirty":identity.git_dirty,"argv":list(argv),"shell_rendering":shlex.join(argv),"working_directory":str(work),"environment_overrides":env_visible,"redacted_environment_keys":list(spec.redacted_environment_keys),"started_utc":started,"stdout_log":str(stdout_path),"stderr_log":str(stderr_path),"combined_log":str(combined_path),"declared_inputs":_serialize(spec.declared_inputs),"declared_outputs":_serialize(spec.declared_outputs),"input_file_identities":[asdict(file_identity(x.path)) for x in spec.declared_inputs if Path(x.path).is_file()],"expected_exit_codes":list(spec.expected_exit_codes),"timeout_seconds":spec.timeout_seconds,"execution_mode":spec.execution_mode.value,"container_identity":container_identity,"dry_run":dry_run,"overwrite":spec.overwrite}
    def finish(status,exit_code=None,failure=None,outputs=()):
        record={**base,"completed_utc":utc_now(),"duration_seconds":round(time.monotonic()-begin,6),"exit_code":exit_code,"output_file_identities":list(outputs),"status":status.value,"failure":failure}; atomic_write_json(record_path,record); return record
    if dry_run:
        message="DRY RUN: command was not executed\n"; stdout_path.write_text(message); stderr_path.write_text(""); combined_path.write_text(message)
        return finish(ExecutionStatus.DRY_RUN)
    atomic_write_json(record_path,{**base,"completed_utc":None,"duration_seconds":None,"exit_code":None,"output_file_identities":[],"status":"RUNNING","failure":None})
    env=os.environ.copy(); env.update(spec.environment_overrides)
    with stdout_path.open("wb") as out, stderr_path.open("wb") as err:
        process=subprocess.Popen(list(argv),cwd=work,env=env,stdout=out,stderr=err,start_new_session=(os.name!="nt"))
        try:
            while True:
                if cancellation and cancellation.cancelled: raise KeyboardInterrupt
                try: exit_code=process.wait(timeout=.05); break
                except subprocess.TimeoutExpired:
                    if spec.timeout_seconds is not None and time.monotonic()-begin >= spec.timeout_seconds: raise TimeoutError
        except (TimeoutError,KeyboardInterrupt) as exc:
            if os.name!="nt": os.killpg(process.pid,signal.SIGTERM)
            else: process.terminate()
            try: process.wait(timeout=termination_grace_seconds)
            except subprocess.TimeoutExpired:
                if os.name!="nt": os.killpg(process.pid,signal.SIGKILL)
                else: process.kill()
                process.wait()
            status=ExecutionStatus.TIMED_OUT if isinstance(exc,TimeoutError) else ExecutionStatus.CANCELLED
            failure={"error_type":"CommandTimeoutError" if status is ExecutionStatus.TIMED_OUT else "Cancelled","message":"command timed out" if status is ExecutionStatus.TIMED_OUT else "command cancelled","exit_code":None,"signal":"SIGTERM","timed_out":status is ExecutionStatus.TIMED_OUT,"missing_outputs":[],"diagnostic_log_tail":""}
            finish(status,None,failure)
            if status is ExecutionStatus.TIMED_OUT: raise CommandTimeoutError(failure["message"],command_id=spec.command_id)
            raise CommandExecutionError(failure["message"],command_id=spec.command_id)
    with combined_path.open("wb") as combined:
        for label,path in ((b"[stdout]\n",stdout_path),(b"[stderr]\n",stderr_path)): combined.write(label); combined.write(path.read_bytes())
    tail=stderr_path.read_bytes()[-8192:].decode("utf-8","replace")
    if exit_code not in spec.expected_exit_codes:
        failure={"error_type":"CommandExecutionError","message":f"unexpected exit code {exit_code}","exit_code":exit_code,"signal":None,"timed_out":False,"missing_outputs":[],"diagnostic_log_tail":tail}
        finish(ExecutionStatus.FAILED,exit_code,failure); raise CommandExecutionError(failure["message"],exit_code=exit_code)
    missing=[o.path for o in spec.declared_outputs if o.required and (not Path(o.path).is_file() or (not o.allow_empty and Path(o.path).stat().st_size==0))]
    if missing:
        failure={"error_type":"ExpectedOutputMissingError","message":"required output missing or empty","exit_code":exit_code,"signal":None,"timed_out":False,"missing_outputs":missing,"diagnostic_log_tail":tail}
        finish(ExecutionStatus.OUTPUT_VALIDATION_FAILED,exit_code,failure); raise ExpectedOutputMissingError(failure["message"],missing_outputs=missing)
    outputs=[asdict(file_identity(o.path)) for o in spec.declared_outputs if Path(o.path).is_file()]
    return finish(ExecutionStatus.SUCCEEDED,exit_code,None,outputs)
