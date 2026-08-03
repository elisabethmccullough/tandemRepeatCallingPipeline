"""Canonical stage registry, lifecycle vocabulary, and conservative resume checks."""
from dataclasses import asdict, dataclass
from enum import Enum
from pathlib import Path
from typing import Mapping
from .provenance import file_identity

class StageStatus(str,Enum):
    NOT_STARTED="NOT_STARTED"; PLANNED="PLANNED"; RUNNING="RUNNING"; SUCCEEDED="SUCCEEDED"; FAILED="FAILED"; SKIPPED="SKIPPED"; DRY_RUN="DRY_RUN"; INVALIDATED="INVALIDATED"
class ResumeReason(str,Enum):
    RESUME_ALLOWED="RESUME_ALLOWED"; NO_PRIOR_RECORD="NO_PRIOR_RECORD"; PRIOR_STAGE_FAILED="PRIOR_STAGE_FAILED"; INPUT_CHANGED="INPUT_CHANGED"; OUTPUT_MISSING="OUTPUT_MISSING"; OUTPUT_CHANGED="OUTPUT_CHANGED"; CONFIGURATION_CHANGED="CONFIGURATION_CHANGED"; TOOL_CHANGED="TOOL_CHANGED"; TOOL_VERSION_CHANGED="TOOL_VERSION_CHANGED"; EXECUTION_MODE_CHANGED="EXECUTION_MODE_CHANGED"; CONTAINER_CHANGED="CONTAINER_CHANGED"; UNSUPPORTED_RECORD_VERSION="UNSUPPORTED_RECORD_VERSION"
@dataclass(frozen=True)
class StageDefinition:
    stage_id:str; order:int; display_name:str; description:str; required_input_roles:tuple[str,...]=(); expected_output_roles:tuple[str,...]=(); required_tools:tuple[str,...]=(); optional_tools:tuple[str,...]=(); supports_dry_run:bool=True; skippable:bool=True; resume_supported:bool=True

_IDS=("00_validate_inputs","01_prepare_bam","02_align_assembly","03_run_vamos_read","04_run_vamos_contig","05_run_straglr","06_prepare_tandem_genotypes","07_run_tandem_genotypes","08_normalize_outputs","09_build_case_package","10_validate_case_package")
STAGES=tuple(StageDefinition(x,i,x.replace("_"," ").title(),f"Pipeline scaffold stage {x}.") for i,x in enumerate(_IDS))
if len({x.stage_id for x in STAGES})!=len(STAGES) or len({x.order for x in STAGES})!=len(STAGES): raise RuntimeError("invalid stage registry")

def select_stages(start:str|None=None,stop:str|None=None):
    index={x.stage_id:i for i,x in enumerate(STAGES)}; start=start or STAGES[0].stage_id; stop=stop or STAGES[-1].stage_id
    if start not in index: raise ValueError(f"unknown start stage: {start}")
    if stop not in index: raise ValueError(f"unknown stop stage: {stop}")
    if index[start]>index[stop]: raise ValueError("start stage follows stop stage")
    return STAGES[index[start]:index[stop]+1]

def resume_eligibility(prior:Mapping|None,configuration_digest:str,input_identities:list[Mapping],output_identities:list[Mapping],tool_identities:list[Mapping]):
    if not prior:return ResumeReason.NO_PRIOR_RECORD
    if prior.get("record_schema_version")!="1.0":return ResumeReason.UNSUPPORTED_RECORD_VERSION
    if prior.get("status")!="SUCCEEDED":return ResumeReason.PRIOR_STAGE_FAILED
    if prior.get("configuration_digest")!=configuration_digest:return ResumeReason.CONFIGURATION_CHANGED
    for expected in prior.get("input_file_identities",[]):
        p=Path(expected["path"])
        if not p.is_file() or asdict(file_identity(p))!=expected:return ResumeReason.INPUT_CHANGED
    for expected in prior.get("output_file_identities",[]):
        p=Path(expected["path"])
        if not p.is_file():return ResumeReason.OUTPUT_MISSING
        if asdict(file_identity(p))!=expected:return ResumeReason.OUTPUT_CHANGED
    old=prior.get("tool_identities",[])
    if [{k:x.get(k) for k in ("tool_id","resolved_executable")} for x in old] != [{k:x.get(k) for k in ("tool_id","resolved_executable")} for x in tool_identities]: return ResumeReason.TOOL_CHANGED
    if [x.get("detected_version") for x in old]!=[x.get("detected_version") for x in tool_identities]:return ResumeReason.TOOL_VERSION_CHANGED
    if [x.get("execution_mode") for x in old]!=[x.get("execution_mode") for x in tool_identities]:return ResumeReason.EXECUTION_MODE_CHANGED
    if [x.get("container_digest") for x in old]!=[x.get("container_digest") for x in tool_identities]:return ResumeReason.CONTAINER_CHANGED
    return ResumeReason.RESUME_ALLOWED
