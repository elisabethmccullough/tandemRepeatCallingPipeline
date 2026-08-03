"""Explicit VAMOS 2.x adapter, native registration, and conservative JSONL parser.

The command spelling is deliberately isolated here and must be confirmed against
the laboratory installation before real use. Unknown versions are never guessed.
"""
from __future__ import annotations
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import re
from typing import Any
from uuid import uuid4

from ..caller_evidence import (ASSEMBLED_CONTIG, DIRECT_SEQUENCE_ASSOCIATION,
    RAW_READS, UNASSIGNED, CallerRunStatus)
from ..caller_outputs import NativeCallerOutput
from ..config import load_locus_config
from ..execution import CommandSpec, InputDeclaration, OutputDeclaration, execute
from ..fasta import read_fasta
from ..provenance import atomic_write_json, file_identity, sha256_file, utc_now
from ..tools import Tool, ToolId, ToolStatus, detect_version, resolve_tool

@dataclass(frozen=True)
class VamosCapabilities:
    supports_bam_input: bool = True
    supports_fasta_or_contig_input: bool = True
    supports_reference_argument: bool = True
    supports_catalog_argument: bool = True
    supports_output_prefix: bool = True
    native_output_types: tuple[str, ...] = ("jsonl", "summary.txt")

@dataclass(frozen=True)
class VamosCommandPlan:
    command_id: str
    argv: tuple[str, ...]
    native_outputs: tuple[str, ...]

class UnsupportedVamosVersion(ValueError): pass
class UnsupportedVamosFormat(ValueError): pass

@dataclass(frozen=True)
class VamosAdapter:
    name: str = "vamos-2.x-provisional"
    capabilities: VamosCapabilities = VamosCapabilities()

    def plan(self, executable: str, *, mode: str, input_path: Path, reference: Path,
             catalog: Path, output_prefix: Path, threads: int, additional_arguments=()) -> VamosCommandPlan:
        mode_flag = "--bam" if mode == "read" else "--fasta"
        argv = (executable, mode_flag, str(input_path), "--reference", str(reference),
                "--catalog", str(catalog), "--output-prefix", str(output_prefix),
                "--threads", str(threads), *tuple(additional_arguments))
        return VamosCommandPlan(f"vamos-{mode}-{output_prefix.parent.name}", argv,
            (str(output_prefix.with_suffix(".jsonl")), str(output_prefix.with_suffix(".summary.txt"))))

def select_adapter(version: str | None) -> VamosAdapter:
    if version and re.fullmatch(r"2(?:\.[0-9A-Za-z_-]+)+", version): return VamosAdapter()
    raise UnsupportedVamosVersion("VAMOS command syntax is supported only for explicitly recognized 2.x versions")

def detect_vamos(executable: str, config_directory: Path, required: bool) -> Tool:
    tool = resolve_tool(Tool(ToolId.VAMOS, "VAMOS", executable, required), config_directory)
    return detect_version(tool, pattern=r"(?i)(?:vamos[^0-9]*)?v?([0-9]+(?:\.[0-9A-Za-z_-]+)+)") if tool.resolved_executable else tool

def parse_native_jsonl(path: Path) -> list[dict[str, Any]]:
    records=[]
    try:
        for line in path.read_text(encoding="utf-8").splitlines():
            if not line.strip(): continue
            value=json.loads(line)
            if not isinstance(value, dict) or "record_id" not in value or "allele_id" not in value: raise UnsupportedVamosFormat("required record_id/allele_id fields absent")
            records.append(value)
    except (UnicodeError, json.JSONDecodeError) as exc: raise UnsupportedVamosFormat(str(exc)) from exc
    if not records: raise UnsupportedVamosFormat("no supported VAMOS JSONL records")
    return records

def _normalized(raw: dict[str, Any], *, config: dict[str, Any], version: str | None,
                source: NativeCallerOutput, analysis_source: str, sequence_id: str | None) -> dict[str, Any]:
    locus=load_locus_config(config["locus_config"])["locus"]
    def native(name): return raw.get(name)
    return {"record_schema_version":"1.0", "record_id":f"vamos-{analysis_source.lower()}-{raw['record_id']}-{raw['allele_id']}",
      "case_id":config["run"]["case_id"], "subject_id":config["run"]["subject_id"], "sample_id":config["run"]["sample_id"], "locus_id":config["run"]["locus_id"],
      "caller":"VAMOS", "caller_version":version, "analysis_source":analysis_source,
      "source_file_id":source.file_id, "source_file_sha256":source.sha256,
      "native_record_identifier":str(raw["record_id"]), "native_allele_identifier":str(raw["allele_id"]),
      "associated_sequence_id":sequence_id, "assignment_state":DIRECT_SEQUENCE_ASSOCIATION if sequence_id else UNASSIGNED,
      "reference_build":native("reference_build") or locus.get("reference_build"), "chromosome":native("chromosome"), "start":native("start"), "end":native("end"),
      "coordinate_convention":native("coordinate_convention"), "reported_motif":native("motif"), "reported_motif_chain":native("motif_chain"),
      "reported_repeat_count":native("repeat_count"), "reported_repeat_length_bp":native("repeat_length_bp"), "supporting_reads":native("supporting_reads"),
      "total_spanning_reads":native("total_spanning_reads"), "quality_state":native("quality_state"), "raw_fields":raw, "normalization_warnings":[]}

def _catalog(config: dict[str, Any], mode: str) -> Path | None:
    locus_path=Path(config["locus_config"]); resources=load_locus_config(locus_path)["caller_resources"]["vamos"]
    value=resources.get(f"{mode}_catalog", resources.get("repeat_catalog"))
    if value is None: return None
    candidate=Path(value).expanduser(); return candidate if candidate.is_absolute() else (locus_path.parent/candidate).resolve()

def _base_metadata(config, stage_id, status, source, sequence_id=None, source_record_id=None, tool=None, catalog=None, warnings=(), failure=None):
    now=utc_now()
    return {"record_schema_version":"1.0", "run_id":str(uuid4()), "stage_id":stage_id, "caller":"VAMOS",
      "caller_version":tool.detected_version if tool else None, "raw_version_output":tool.raw_version_output if tool else None,
      "analysis_source":source, "associated_sequence_id":sequence_id, "source_fasta_record_id":source_record_id,
      "assignment_state":DIRECT_SEQUENCE_ASSOCIATION if sequence_id else UNASSIGNED, "status":status.value,
      "command_record_paths":[], "input_file_ids":[], "native_output_file_ids":[],
      "catalog_identity":asdict(file_identity(catalog)) if catalog and catalog.is_file() else None,
      "started_utc":now, "completed_utc":now, "warnings":list(warnings), "failure":failure}

def _run_one(config, stage_id, mode, input_path, out, tool, catalog, *, sequence_id=None, source_record_id=None, overwrite=False):
    source=RAW_READS if mode=="read" else ASSEMBLED_CONTIG
    adapter=select_adapter(tool.detected_version); native=out/"native"; native.mkdir(parents=True,exist_ok=True)
    prefix=native/"vamos"; plan=adapter.plan(tool.resolved_executable, mode=mode, input_path=input_path,
        reference=Path(config["inputs"]["reference_fasta"]), catalog=catalog, output_prefix=prefix,
        threads=config["execution"]["threads"], additional_arguments=config.get("caller_resources",{}).get("vamos",{}).get("additional_arguments",[]))
    recdir=out/"execution-records"; recdir.mkdir(parents=True,exist_ok=True); record_path=recdir/f"{plan.command_id}.json"
    spec=CommandSpec(plan.command_id,stage_id,"VAMOS",plan.argv,str(out.resolve()),declared_inputs=(InputDeclaration("caller-input",str(input_path)),InputDeclaration("reference",config["inputs"]["reference_fasta"]),InputDeclaration("catalog",str(catalog))),declared_outputs=tuple(OutputDeclaration(f"native-{i}",p,required=True,media_type=None) for i,p in enumerate(plan.native_outputs)),overwrite=overwrite)
    metadata=_base_metadata(config,stage_id,CallerRunStatus.RUNNING,source,sequence_id,source_record_id,tool,catalog)
    metadata["input_file_ids"]=[f"input-{sha256_file(p)[:16]}" for p in (input_path,Path(config["inputs"]["reference_fasta"]),catalog)]
    try:
        execute(spec,tool,record_path,out/"logs",dry_run=False)
        outputs=[NativeCallerOutput.from_path(p,file_id=f"{plan.command_id}-native-{i+1}",caller_version=tool.detected_version,analysis_source=source,producer_command_id=plan.command_id) for i,p in enumerate(plan.native_outputs)]
        normalized=[]; warnings=[]
        try:
            for item in parse_native_jsonl(Path(plan.native_outputs[0])): normalized.append(_normalized(item,config=config,version=tool.detected_version,source=outputs[0],analysis_source=source,sequence_id=sequence_id))
            status=CallerRunStatus.SUCCEEDED
        except UnsupportedVamosFormat as exc:
            status=CallerRunStatus.UNSUPPORTED_FORMAT; warnings=[str(exc)]
        metadata.update(status=status.value,command_record_paths=[str(record_path)],native_output_file_ids=[x.file_id for x in outputs],completed_utc=utc_now(),warnings=warnings)
        atomic_write_json(out/("vamos-read.outputs.json" if mode=="read" else "outputs.json"),{"record_schema_version":"1.0","outputs":[x.to_dict() for x in outputs]})
        atomic_write_json(out/("vamos-read.normalized.json" if mode=="read" else "normalized.json"),{"record_schema_version":"1.0","evidence_state":"AVAILABLE" if normalized else "UNSUPPORTED_FORMAT","records":normalized,"normalization_warnings":warnings})
    except Exception as exc:
        metadata.update(status=CallerRunStatus.FAILED.value,completed_utc=utc_now(),failure={"error_type":type(exc).__name__,"message":str(exc)})
        atomic_write_json(out/("vamos-read.run.json" if mode=="read" else "run.json"),metadata); raise
    atomic_write_json(out/("vamos-read.run.json" if mode=="read" else "run.json"),metadata)
    return metadata, outputs, normalized

def run_vamos_stage(stage_id: str, config: dict[str, Any], root: Path, config_directory: Path, *, overwrite=False, dry_run=False):
    settings=config.get("tools",{}).get("vamos",{}); required=settings.get("required",False)
    tool=detect_vamos(settings.get("executable","vamos"),config_directory,required)
    mode="read" if stage_id=="03_run_vamos_read" else "contig"; source=RAW_READS if mode=="read" else ASSEMBLED_CONTIG
    out=root/("04_vamos_read" if mode=="read" else "05_vamos_contig"); out.mkdir(parents=True,exist_ok=True)
    catalog=_catalog(config,mode)
    if dry_run:
        warning = "Dry run: no VAMOS command or native biological output was created."
        if catalog is None or not catalog.is_file(): warning += " Required catalog is missing."
        if not tool.resolved_executable: warning += " Configured VAMOS executable is unavailable."
        metadata=_base_metadata(config,stage_id,CallerRunStatus.DRY_RUN,source,tool=tool,catalog=catalog,warnings=(warning,))
        atomic_write_json(out/("vamos-read.run.json" if mode=="read" else "stage-summary.json"), metadata if mode=="read" else {"record_schema_version":"1.0","stage_id":stage_id,"runs":[metadata]})
        atomic_write_json(out/("vamos-read.outputs.json" if mode=="read" else "stage-outputs.json"),{"record_schema_version":"1.0","outputs":[]})
        atomic_write_json(out/("vamos-read.normalized.json" if mode=="read" else "stage-normalized.json"),{"record_schema_version":"1.0","evidence_state":"NOT_COMPUTED","records":[],"normalization_warnings":[warning]})
        return
    terminal=None
    if not tool.resolved_executable: terminal=CallerRunStatus.TOOL_MISSING
    elif tool.status is ToolStatus.VERSION_UNDETERMINED: terminal=CallerRunStatus.UNSUPPORTED_VERSION
    else:
        try: select_adapter(tool.detected_version)
        except UnsupportedVamosVersion: terminal=CallerRunStatus.UNSUPPORTED_VERSION
    if catalog is None or not catalog.is_file(): terminal=CallerRunStatus.INPUT_MISSING
    if terminal:
        metadata=_base_metadata(config,stage_id,terminal,source,tool=tool,catalog=catalog,warnings=("VAMOS evidence was not computed.",))
        target=out/("vamos-read.run.json" if mode=="read" else "stage-summary.json"); atomic_write_json(target,metadata)
        atomic_write_json(out/("vamos-read.outputs.json" if mode=="read" else "stage-outputs.json"),{"record_schema_version":"1.0","outputs":[]})
        atomic_write_json(out/("vamos-read.normalized.json" if mode=="read" else "stage-normalized.json"),{"record_schema_version":"1.0","evidence_state":"INPUT_MISSING" if terminal is CallerRunStatus.INPUT_MISSING else "NOT_COMPUTED","records":[],"normalization_warnings":[terminal.value]})
        if required: raise RuntimeError(f"required VAMOS unavailable: {terminal.value}")
        return
    if mode=="read":
        _run_one(config,stage_id,mode,root/"02_prepared_bam"/"prepared.mini.bam",out,tool,catalog,overwrite=overwrite); return
    fasta={r.identifier:r.sequence for r in read_fasta(root/"03_assembly_alignment"/"patient-sequences.fasta")}
    runs=[]; all_outputs=[]; all_records=[]
    for selected in config["inputs"]["assembly_records"]:
        sequence_id=selected["sequence_id"]; seqout=out/sequence_id; input_dir=seqout/"input"; input_dir.mkdir(parents=True,exist_ok=True)
        sequence=fasta[sequence_id]; single=input_dir/f"{sequence_id}.fasta"
        content=f">{sequence_id}\n{sequence}\n"; single.write_text(content,encoding="utf-8")
        metadata,outputs,records=_run_one(config,stage_id,mode,single,seqout,tool,catalog,sequence_id=sequence_id,source_record_id=selected["record_id"],overwrite=overwrite)
        metadata["sequence_sha256"]=sha256_file(single); atomic_write_json(seqout/"run.json",metadata)
        runs.append(metadata); all_outputs.extend(x.to_dict() for x in outputs); all_records.extend(records)
    atomic_write_json(out/"stage-summary.json",{"record_schema_version":"1.0","stage_id":stage_id,"runs":runs})
    atomic_write_json(out/"stage-outputs.json",{"record_schema_version":"1.0","outputs":all_outputs})
    atomic_write_json(out/"stage-normalized.json",{"record_schema_version":"1.0","evidence_state":"AVAILABLE","records":all_records,"normalization_warnings":[]})
