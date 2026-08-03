"""Stable tool identities, resolution, version probing, and container commands."""
from __future__ import annotations
from dataclasses import asdict, dataclass, replace
from enum import Enum
from pathlib import Path
import os, re, shutil, subprocess
from .provenance import sha256_file

class ToolId(str, Enum):
    SAMTOOLS="SAMTOOLS"; MINIMAP2="MINIMAP2"; VAMOS="VAMOS"; STRAGLR="STRAGLR"; LASTDB="LASTDB"; LASTAL="LASTAL"; TANDEM_GENOTYPES="TANDEM_GENOTYPES"; PYTHON="PYTHON"; PIPELINE="PIPELINE"; APPTAINER="APPTAINER"
class ToolStatus(str, Enum):
    AVAILABLE="AVAILABLE"; MISSING_OPTIONAL="MISSING_OPTIONAL"; MISSING_REQUIRED="MISSING_REQUIRED"; VERSION_UNDETERMINED="VERSION_UNDETERMINED"; UNSUPPORTED_VERSION="UNSUPPORTED_VERSION"; NOT_CHECKED="NOT_CHECKED"
class ExecutionMode(str, Enum): NATIVE="NATIVE"; APPTAINER="APPTAINER"

@dataclass(frozen=True)
class Tool:
    tool_id: ToolId; display_name: str; configured_executable: str; required: bool=False
    resolved_executable: str|None=None; detected_version: str|None=None; raw_version_output: str|None=None
    execution_mode: ExecutionMode=ExecutionMode.NATIVE; status: ToolStatus=ToolStatus.NOT_CHECKED; status_message: str|None=None
    def to_dict(self): return {k:(v.value if isinstance(v, Enum) else v) for k,v in asdict(self).items()}

def resolve_tool(tool: Tool, config_directory: str|Path, path: str|None=None) -> Tool:
    configured = Path(tool.configured_executable).expanduser()
    configured_text = tool.configured_executable
    explicit = configured.is_absolute() or configured.parent != Path(".") or os.sep in configured_text or bool(os.altsep and os.altsep in configured_text)
    candidate = str((Path(config_directory)/configured).resolve()) if explicit and not configured.is_absolute() else str(configured)
    resolved = shutil.which(candidate, path=path) if not explicit else (str(Path(candidate).resolve()) if Path(candidate).is_file() else None)
    if resolved is None:
        status = ToolStatus.MISSING_REQUIRED if tool.required else ToolStatus.MISSING_OPTIONAL
        return replace(tool, status=status, status_message="configured executable was not found")
    return replace(tool, resolved_executable=resolved, status=ToolStatus.AVAILABLE)

def detect_version(tool: Tool, arguments: tuple[str,...]=( "--version",), pattern: str=r"(?i)(?:version\s*)?v?([0-9]+(?:\.[0-9A-Za-z_-]+)+)", accepted_exit_codes: tuple[int,...]=(0,1)) -> Tool:
    if not tool.resolved_executable: return tool
    result=subprocess.run([tool.resolved_executable,*arguments],capture_output=True,text=True,errors="replace")
    raw=result.stdout+result.stderr; match=re.search(pattern, raw)
    if result.returncode in accepted_exit_codes and match: return replace(tool,detected_version=match.group(1),raw_version_output=raw,status=ToolStatus.AVAILABLE)
    return replace(tool,raw_version_output=raw,status=ToolStatus.VERSION_UNDETERMINED,status_message="version output could not be parsed")

@dataclass(frozen=True)
class ContainerMetadata:
    apptainer_executable: str; container_image: str; bind_mounts: tuple[str,...]=(); internal_executable: str=""; container_working_directory: str|None=None
    @property
    def container_sha256(self): return sha256_file(self.container_image) if Path(self.container_image).is_file() else None
    def host_argv(self, internal_argv: tuple[str,...], host_working_directory: str) -> tuple[str,...]:
        args=[self.apptainer_executable,"exec"]
        for bind in self.bind_mounts: args += ["--bind",bind]
        if self.container_working_directory: args += ["--pwd",self.container_working_directory]
        return tuple(args+[self.container_image,self.internal_executable,*internal_argv])
