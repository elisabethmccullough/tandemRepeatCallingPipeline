"""Conservative, development-gated STRaglr execution and normalization.

No real STRaglr release is asserted as verified here.  The sole adapter is for
synthetic integration tests and must be explicitly enabled by configuration.
"""
from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
from enum import Enum
import os
from pathlib import Path
import re
import shutil
from typing import Any
from uuid import uuid4

from ..caller_evidence import CallerRunStatus, RAW_READS, UNASSIGNED
from ..caller_outputs import NativeCallerOutput
from ..config import load_locus_config
from ..execution import CommandSpec, InputDeclaration, OutputDeclaration, execute
from ..provenance import atomic_write_json, file_identity, sha256_file, utc_now
from ..tools import Tool, ToolId, ToolStatus, detect_version, resolve_tool


class StraglrVersionClassification(str, Enum):
    VERIFIED_SUPPORTED = "VERIFIED_SUPPORTED"
    PROVISIONAL_DEVELOPMENT = "PROVISIONAL_DEVELOPMENT"
    UNSUPPORTED = "UNSUPPORTED"
    UNDETERMINED = "UNDETERMINED"


@dataclass(frozen=True)
class StraglrVersionRange:
    major: int
    classification: StraglrVersionClassification


@dataclass(frozen=True)
class StraglrCapabilities:
    native_coordinate_convention: str = "zero_based_half_open"
    native_output_types: tuple[str, ...] = ("tsv",)
    supports_threads: bool = True


@dataclass(frozen=True)
class StraglrCommandPlan:
    command_id: str
    argv: tuple[str, ...]
    native_outputs: tuple[str, ...]


@dataclass(frozen=True)
class StraglrAdapter:
    name: str = "straglr-1.x-development-provisional"
    version_range: StraglrVersionRange = StraglrVersionRange(1, StraglrVersionClassification.PROVISIONAL_DEVELOPMENT)
    capabilities: StraglrCapabilities = StraglrCapabilities()

    def plan(self, executable: str, *, bam: Path, reference: Path, catalog: Path,
             output_prefix: Path, threads: int, additional_arguments=()) -> StraglrCommandPlan:
        output = output_prefix.with_suffix(".tsv")
        # This spelling is intentionally provisional, isolated, and not a claim
        # about an installed laboratory STRaglr release.
        argv = (executable, str(bam), str(reference), str(catalog), str(output),
                "--threads", str(threads), *tuple(additional_arguments))
        return StraglrCommandPlan("straglr-read", argv, (str(output),))


@dataclass(frozen=True)
class StraglrNativeRecord:
    native_record_identifier: str
    native_allele_identifier: str
    row_number: int
    raw_fields: dict[str, str]


class UnsupportedStraglrVersion(ValueError): pass
class UnsupportedStraglrFormat(ValueError): pass


def classify_version(version: str | None) -> StraglrVersionClassification:
    if version is None:
        return StraglrVersionClassification.UNDETERMINED
    if re.fullmatch(r"1(?:\.[0-9A-Za-z_-]+)+", version):
        return StraglrVersionClassification.PROVISIONAL_DEVELOPMENT
    return StraglrVersionClassification.UNSUPPORTED


def select_adapter(version: str | None, *, allow_provisional: bool = False) -> StraglrAdapter:
    classification = classify_version(version)
    if classification is StraglrVersionClassification.PROVISIONAL_DEVELOPMENT and allow_provisional:
        return StraglrAdapter()
    if classification is StraglrVersionClassification.PROVISIONAL_DEVELOPMENT:
        raise UnsupportedStraglrVersion("STRaglr has no laboratory-verified adapter; provisional adapter is disabled")
    if classification is StraglrVersionClassification.UNDETERMINED:
        raise UnsupportedStraglrVersion("STRaglr version could not be determined; command syntax will not be guessed")
    raise UnsupportedStraglrVersion(f"unsupported STRaglr version: {version}")


def detect_straglr(executable: str, config_directory: Path, required: bool) -> Tool:
    tool = resolve_tool(Tool(ToolId.STRAGLR, "STRaglr", executable, required), config_directory)
    if not tool.resolved_executable:
        return tool
    return detect_version(tool, pattern=r"(?i)(?:straglr[^0-9]*)?v?([0-9]+(?:\.[0-9A-Za-z_-]+)+)")


REQUIRED_COLUMNS = ("locus_id", "chromosome", "start", "end", "allele_id")


def parse_native_tsv(path: str | Path) -> list[StraglrNativeRecord]:
    """Parse the frozen development TSV layout without coercing any values."""
    source = Path(path)
    try:
        with source.open(encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            if reader.fieldnames is None or any(name not in reader.fieldnames for name in REQUIRED_COLUMNS):
                raise UnsupportedStraglrFormat("unsupported STRaglr TSV columns")
            if None in reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
                raise UnsupportedStraglrFormat("ambiguous or duplicate STRaglr TSV columns")
            records = []
            for row_number, row in enumerate(reader, 1):
                if None in row or any(value is None for value in row.values()):
                    raise UnsupportedStraglrFormat(f"malformed STRaglr row {row_number}: column count differs from header")
                allele = row["allele_id"]
                if not allele:
                    raise UnsupportedStraglrFormat(f"malformed STRaglr row {row_number}: allele_id is empty")
                native_id = row.get("record_id") or f"row-{row_number}"
                records.append(StraglrNativeRecord(native_id, allele, row_number, dict(row)))
    except UnicodeError as exc:
        raise UnsupportedStraglrFormat(str(exc)) from exc
    if not records:
        raise UnsupportedStraglrFormat("no supported STRaglr records")
    return records


def _integer(raw: str | None, field: str, warnings: list[str]) -> int | None:
    if raw in (None, "", ".", "NA"):
        return None
    try:
        value = int(raw)
        if value < 0: raise ValueError
        return value
    except ValueError:
        warnings.append(f"{field} is not a non-negative integer; normalized value is null")
        return None


def normalize_record(native: StraglrNativeRecord, *, config: dict[str, Any], caller_version: str | None,
                     source: NativeCallerOutput, coordinate_verified: bool = True) -> dict[str, Any]:
    raw = native.raw_fields
    warnings: list[str] = []
    start = _integer(raw.get("start"), "start", warnings) if coordinate_verified else None
    end = _integer(raw.get("end"), "end", warnings) if coordinate_verified else None
    if not coordinate_verified:
        warnings.append("native coordinate convention is unverified; normalized coordinates are null")
    if start is not None and end is not None and start > end:
        warnings.append("native end precedes start; normalized coordinates are null")
        start = end = None
    locus = load_locus_config(config["locus_config"])["locus"]
    return {
        "record_schema_version":"1.0", "record_id":f"straglr-{native.native_record_identifier}-{native.native_allele_identifier}",
        "case_id":config["run"]["case_id"], "subject_id":config["run"]["subject_id"], "sample_id":config["run"]["sample_id"],
        "locus_id":config["run"]["locus_id"], "caller":"STRAGLR", "caller_version":caller_version,
        "analysis_source":RAW_READS, "source_file_id":source.file_id, "source_file_sha256":source.sha256,
        "native_record_identifier":native.native_record_identifier, "native_allele_identifier":native.native_allele_identifier,
        "associated_sequence_id":None, "assignment_state":UNASSIGNED, "reference_build":raw.get("reference_build") or locus["reference_build"],
        "chromosome":raw.get("chromosome") or None, "start":start, "end":end,
        "coordinate_convention":"zero_based_half_open", "reported_motif":raw.get("repeat_unit") or None,
        "reported_motif_chain":None, "reported_repeat_count":raw.get("repeat_count") or None,
        "reported_repeat_length_bp":raw.get("repeat_size") or None,
        "supporting_reads":_integer(raw.get("supporting_reads"), "supporting_reads", warnings),
        "total_spanning_reads":_integer(raw.get("total_spanning_reads"), "total_spanning_reads", warnings),
        "quality_state":"AVAILABLE" if not warnings else "AMBIGUOUS", "raw_fields":raw,
        "normalization_warnings":warnings,
    }


def _catalog(config: dict[str, Any]) -> Path | None:
    locus_path = Path(config["locus_config"])
    value = _resources(config).get("repeat_catalog")
    if value is None: return None
    candidate = Path(value).expanduser()
    return candidate if candidate.is_absolute() else (locus_path.parent / candidate).resolve()


def _resources(config: dict[str, Any]) -> dict[str, Any]:
    """Return the single canonical STRaglr caller-resource configuration."""
    return load_locus_config(config["locus_config"])["caller_resources"]["straglr"]


def _metadata(config, status, tool, catalog, *, warnings=(), failure=None):
    now = utc_now()
    return {"record_schema_version":"1.0", "run_id":str(uuid4()), "stage_id":"05_run_straglr",
        "caller":"STRAGLR", "caller_version":tool.detected_version if tool else None,
        "raw_version_output":tool.raw_version_output if tool else None, "analysis_source":RAW_READS,
        "associated_sequence_id":None, "assignment_state":UNASSIGNED, "status":status.value,
        "command_record_paths":[], "input_file_ids":[], "native_output_file_ids":[],
        "catalog_identity":asdict(file_identity(catalog)) if catalog and catalog.is_file() else None,
        "caller_resource_settings":{"allow_provisional_adapter":bool(_resources(config).get("allow_provisional_adapter",False)),
            "additional_arguments":list(_resources(config).get("additional_arguments",[]))},
        "started_utc":now, "completed_utc":now, "warnings":list(warnings), "failure":failure}


def run_straglr_stage(config: dict[str, Any], root: Path, config_directory: Path, *, overwrite=False, dry_run=False):
    out=root/"06_straglr"; out.mkdir(parents=True,exist_ok=True)
    settings=config.get("tools",{}).get("straglr",{}); required=settings.get("required",False)
    resources=_resources(config)
    tool=detect_straglr(settings.get("executable","straglr.py"),config_directory,required)
    catalog=_catalog(config)
    def terminal(status, warning, evidence="NOT_COMPUTED"):
        metadata=_metadata(config,status,tool,catalog,warnings=(warning,))
        atomic_write_json(out/"straglr.run.json",metadata)
        atomic_write_json(out/"straglr.outputs.json",{"record_schema_version":"1.0","outputs":[]})
        atomic_write_json(out/"straglr.normalized.json",{"record_schema_version":"1.0","evidence_state":evidence,"records":[],"normalization_warnings":[warning]})
        if required and not dry_run: raise RuntimeError(f"required STRaglr unavailable: {status.value}")
        return metadata
    if dry_run:
        details=[]
        if not tool.resolved_executable: details.append("configured executable unavailable")
        if catalog is None or not catalog.is_file(): details.append("required catalog missing")
        return terminal(CallerRunStatus.DRY_RUN,"Dry run: no STRaglr command or biological output was created"+(f"; {', '.join(details)}" if details else ""))
    if not tool.resolved_executable: return terminal(CallerRunStatus.TOOL_MISSING,"configured STRaglr executable is unavailable")
    if tool.status is ToolStatus.VERSION_UNDETERMINED: return terminal(CallerRunStatus.UNSUPPORTED_VERSION,"STRaglr version is undetermined")
    try: adapter=select_adapter(tool.detected_version,allow_provisional=resources.get("allow_provisional_adapter",False))
    except UnsupportedStraglrVersion as exc: return terminal(CallerRunStatus.UNSUPPORTED_VERSION,str(exc))
    if catalog is None or not catalog.is_file(): return terminal(CallerRunStatus.INPUT_MISSING,"explicit STRaglr repeat catalog is missing","INPUT_MISSING")
    native=out/"native"
    if native.exists() and not overwrite:
        raise FileExistsError(f"completed STRaglr native output set already exists: {native}")
    work_native=out/f".native-work-{uuid4()}"; work_native.mkdir(parents=True)
    bam=root/"02_prepared_bam"/"prepared.mini.bam"; bai=root/"02_prepared_bam"/"prepared.mini.bam.bai"
    reference=Path(config["inputs"]["reference_fasta"]); reference_index=Path(config["inputs"]["reference_fasta_index"])
    plan=adapter.plan(tool.resolved_executable,bam=bam,reference=reference,catalog=catalog,output_prefix=work_native/"straglr",
        threads=config["execution"]["threads"],additional_arguments=resources.get("additional_arguments",[]))
    record_path=out/"execution-records"/f"{plan.command_id}.json"
    declarations=(InputDeclaration("prepared-bam",str(bam)),InputDeclaration("prepared-bam-index",str(bai)),
        InputDeclaration("reference",str(reference)),InputDeclaration("reference-index",str(reference_index)),InputDeclaration("catalog",str(catalog)))
    spec=CommandSpec(plan.command_id,"05_run_straglr","STRAGLR",plan.argv,str(out.resolve()),declared_inputs=declarations,
        declared_outputs=tuple(OutputDeclaration(f"native-{i}",p,required=True,media_type=None) for i,p in enumerate(plan.native_outputs)),overwrite=overwrite)
    warning="Development-only provisional STRaglr adapter enabled; command syntax and TSV layout are not laboratory verified."
    metadata=_metadata(config,CallerRunStatus.RUNNING,tool,catalog,warnings=(warning,)); metadata["input_file_ids"]=[f"input-{sha256_file(x.path)[:16]}" for x in declarations]
    try:
        execute(spec,tool,record_path,out/"logs")
        produced_names=sorted(path.name for path in work_native.iterdir() if path.is_file())
        if not produced_names:
            raise UnsupportedStraglrFormat("STRaglr execution produced no native files")
        # Publish one completed execution as a directory replacement. Old files
        # never participate in discovery and failed work is never published.
        old_native=out/f".native-old-{uuid4()}"
        if native.exists(): os.replace(native,old_native)
        try:
            os.replace(work_native,native)
        except Exception:
            if old_native.exists(): os.replace(old_native,native)
            raise
        if old_native.exists(): shutil.rmtree(old_native)
        native_paths=[native/name for name in produced_names]
        outputs=[NativeCallerOutput.from_path(path,file_id=f"straglr-native-{i+1}",caller="STRAGLR",caller_version=tool.detected_version,
            analysis_source=RAW_READS,producer_command_id=plan.command_id) for i,path in enumerate(native_paths)]
        records=[]; normalization_warnings=[warning]; status=CallerRunStatus.SUCCEEDED
        try:
            primary=native/Path(plan.native_outputs[0]).name
            primary_output=next(item for item in outputs if Path(item.path)==primary)
            records=[normalize_record(row,config=config,caller_version=tool.detected_version,source=primary_output) for row in parse_native_tsv(primary)]
        except UnsupportedStraglrFormat as exc:
            status=CallerRunStatus.UNSUPPORTED_FORMAT; normalization_warnings.append(str(exc))
        metadata.update(status=status.value,command_record_paths=[str(record_path)],native_output_file_ids=[x.file_id for x in outputs],completed_utc=utc_now(),warnings=normalization_warnings)
        atomic_write_json(out/"straglr.outputs.json",{"record_schema_version":"1.0","outputs":[x.to_dict() for x in outputs]})
        atomic_write_json(out/"straglr.normalized.json",{"record_schema_version":"1.0","evidence_state":"AVAILABLE" if records else "UNSUPPORTED_FORMAT","records":records,"normalization_warnings":normalization_warnings})
    except Exception as exc:
        if work_native.exists(): shutil.rmtree(work_native)
        metadata.update(status=CallerRunStatus.FAILED.value,completed_utc=utc_now(),failure={"error_type":type(exc).__name__,"message":str(exc)})
        atomic_write_json(out/"straglr.run.json",metadata); raise
    atomic_write_json(out/"straglr.run.json",metadata)
    return metadata
