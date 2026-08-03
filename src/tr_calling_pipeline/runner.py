"""Scaffold runner. Caller-specific execution intentionally remains future work."""
from dataclasses import asdict
from pathlib import Path
from .config import load_config
from .provenance import atomic_write_json, canonical_digest, utc_now
from .stages import select_stages

def run(config_path,*,dry_run=False,start_stage=None,stop_stage=None,resume=True,overwrite=False,execution_mode=None):
    config=load_config(config_path); stages=select_stages(start_stage,stop_stage)
    root=Path(config["run"]["output_root"])/f'{config["run"]["sample_id"]}_{config["run"]["locus_id"]}'
    records=root/"00_manifest"/"stages"; records.mkdir(parents=True,exist_ok=True)
    for stage in stages:
        now=utc_now(); record={"record_schema_version":"1.0","stage_id":stage.stage_id,"status":"DRY_RUN" if dry_run else "PLANNED","started_utc":now,"completed_utc":now,"duration_seconds":0.0,"configuration_digest":canonical_digest({"stage_id":stage.stage_id,"run":config["run"],"inputs":config["inputs"],"tools":config.get("tools",{}),"container":config.get("container",{})}),"input_file_identities":[],"output_file_identities":[],"tool_identities":[],"command_record_paths":[],"warnings":["Caller-specific execution is not implemented; this is a scaffold plan."],"failure":None,"resume_eligibility":{"eligible":False,"reason":"NO_PRIOR_RECORD"}}
        atomic_write_json(records/f"{stage.stage_id}.json",record)
    return root
