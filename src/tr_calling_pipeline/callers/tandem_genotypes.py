"""Lossless parser/normalizer for the provisional tandem-genotypes TSV contract."""
from __future__ import annotations
import csv, re, json, os, shutil
from dataclasses import asdict, dataclass
from enum import Enum
from pathlib import Path
from typing import Any
from ..caller_outputs import NativeCallerOutput
from ..tools import Tool, ToolId, ToolStatus, detect_version, resolve_tool
from ..execution import CommandSpec, InputDeclaration, OutputDeclaration, execute
from ..provenance import atomic_write_json, file_identity, utc_now
from ..config import load_locus_config
from uuid import uuid4

class TandemGenotypesVersionClassification(str, Enum):
    VERIFIED_SUPPORTED="VERIFIED_SUPPORTED"; PROVISIONAL_DEVELOPMENT="PROVISIONAL_DEVELOPMENT"; UNSUPPORTED="UNSUPPORTED"; UNDETERMINED="UNDETERMINED"
class UnsupportedTandemGenotypesVersion(ValueError): pass
class UnsupportedTandemGenotypesFormat(ValueError): pass

@dataclass(frozen=True)
class TandemGenotypesCapabilities:
    native_output_format: str="development-tsv-v1"
    coordinate_space: str="UNKNOWN_COORDINATE_SPACE"
@dataclass(frozen=True)
class TandemGenotypesCommandPlan:
    command_id: str; argv: tuple[str,...]; native_outputs: tuple[str,...]
@dataclass(frozen=True)
class TandemGenotypesAdapter:
    name: str="tandem-genotypes-development-provisional"
    capabilities: TandemGenotypesCapabilities=TandemGenotypesCapabilities()
    def plan(self, executable: str, *, alignment: Path, repeat_definition: Path, output: Path, additional_arguments=()):
        return TandemGenotypesCommandPlan("tandem-genotypes", (executable, str(alignment), str(repeat_definition), str(output), *tuple(additional_arguments)), (str(output),))
@dataclass(frozen=True)
class TandemGenotypesNativeRecord:
    native_record_identifier: str; native_allele_identifier: str; row_number: int; raw_fields: dict[str,str]; original_text: str

def classify_version(version: str|None):
    if version is None: return TandemGenotypesVersionClassification.UNDETERMINED
    if re.fullmatch(r"0\.[0-9]+(?:\.[0-9]+)?(?:[-+][0-9A-Za-z.-]+)?",version): return TandemGenotypesVersionClassification.PROVISIONAL_DEVELOPMENT
    return TandemGenotypesVersionClassification.UNSUPPORTED
def select_adapter(version: str|None, *, allow_provisional=False):
    kind=classify_version(version)
    if kind is TandemGenotypesVersionClassification.PROVISIONAL_DEVELOPMENT and allow_provisional: return TandemGenotypesAdapter()
    if kind is TandemGenotypesVersionClassification.UNDETERMINED: raise UnsupportedTandemGenotypesVersion("version undetermined; syntax will not be guessed")
    if kind is TandemGenotypesVersionClassification.PROVISIONAL_DEVELOPMENT: raise UnsupportedTandemGenotypesVersion("provisional adapter is disabled")
    raise UnsupportedTandemGenotypesVersion(f"unsupported tandem-genotypes version: {version}")

REQUIRED_COLUMNS=("record_id","allele_id","sequence_id")
def parse_native_tsv(path: str|Path) -> list[TandemGenotypesNativeRecord]:
    lines=Path(path).read_text(encoding="utf-8").splitlines()
    if not lines: raise UnsupportedTandemGenotypesFormat("native output is empty")
    reader=csv.DictReader(lines,delimiter="\t")
    if reader.fieldnames is None or any(x not in reader.fieldnames for x in REQUIRED_COLUMNS) or len(reader.fieldnames)!=len(set(reader.fieldnames)):
        raise UnsupportedTandemGenotypesFormat("unsupported tandem-genotypes TSV columns")
    out=[]
    for number,row in enumerate(reader,1):
        if None in row or any(v is None for v in row.values()): raise UnsupportedTandemGenotypesFormat(f"malformed row {number}")
        if not all(row[x] for x in REQUIRED_COLUMNS): raise UnsupportedTandemGenotypesFormat(f"missing required value in row {number}")
        out.append(TandemGenotypesNativeRecord(row["record_id"],row["allele_id"],number,dict(row),lines[number]))
    if not out: raise UnsupportedTandemGenotypesFormat("native output has no records")
    return out

def normalize_record(native: TandemGenotypesNativeRecord, *, config: dict[str,Any], caller_version: str|None,
                     source: NativeCallerOutput, associated_sequence_id: str) -> dict[str,Any]:
    if native.raw_fields["sequence_id"] != associated_sequence_id: raise ValueError("native record belongs to a different patient sequence")
    r=native.raw_fields
    warnings=["native coordinate space is unverified; normalized coordinates are null"]
    return {"record_schema_version":"1.0","record_id":f"tandem-genotypes-{native.native_record_identifier}-{native.native_allele_identifier}",
      "case_id":config["run"]["case_id"],"subject_id":config["run"]["subject_id"],"sample_id":config["run"]["sample_id"],"locus_id":config["run"]["locus_id"],
      "caller":"TANDEM_GENOTYPES","caller_version":caller_version,"analysis_source":"ASSEMBLED_CONTIG","source_file_id":source.file_id,"source_file_sha256":source.sha256,
      "native_record_identifier":native.native_record_identifier,"native_allele_identifier":native.native_allele_identifier,"associated_sequence_id":associated_sequence_id,
      "assignment_state":"DIRECT_SEQUENCE_ASSOCIATION","reference_build":r.get("reference_build") or None,"chromosome":r.get("chromosome") or None,
      "start":None,"end":None,"coordinate_convention":None,"coordinate_space":"UNKNOWN_COORDINATE_SPACE","reported_motif":r.get("motif") or None,
      "reported_motif_chain":r.get("motif_chain") or None,"reported_repeat_count":r.get("repeat_count") or None,"reported_repeat_length_bp":r.get("repeat_length_bp") or None,
      "supporting_reads":None,"total_spanning_reads":None,"quality_state":"AMBIGUOUS","raw_fields":r,"normalization_warnings":warnings}

def _resource(config):
    locus=Path(config["locus_config"]); value=load_locus_config(locus)["caller_resources"]["tandem_genotypes"]["repeat_definition"]
    if value is None:return None
    p=Path(value).expanduser(); return p if p.is_absolute() else (locus.parent/p).resolve()

def run_tandem_genotypes_stage(config, root: Path, config_directory: Path, *, overwrite=False, dry_run=False):
    out=root/"08_tandem_genotypes"; out.mkdir(parents=True,exist_ok=True)
    tc=config.get("tools",{}).get("tandem_genotypes",{}); required=tc.get("required",False)
    tool=resolve_tool(Tool(ToolId.TANDEM_GENOTYPES,"tandem-genotypes",tc.get("executable","tandem-genotypes"),required),config_directory)
    if tool.resolved_executable: tool=detect_version(tool,pattern=r"(?i)(?:tandem-genotypes[^0-9]*)?v?([0-9]+(?:\.[0-9]+)+)")
    resources=load_locus_config(config["locus_config"])["caller_resources"]["tandem_genotypes"]; repeat=_resource(config)
    prep_path=root/"07_tandem_genotypes_preparation"/"alignment-metadata.json"; preparations=json.loads(prep_path.read_text())["records"] if prep_path.is_file() else []
    terminal=None
    if dry_run: terminal="DRY_RUN"
    elif not tool.resolved_executable: terminal="TOOL_MISSING"
    elif tool.status is ToolStatus.VERSION_UNDETERMINED: terminal="UNSUPPORTED_VERSION"
    elif repeat is None or not repeat.is_file(): terminal="INPUT_MISSING"
    try: adapter=select_adapter(tool.detected_version,allow_provisional=resources["allow_provisional_tandem_genotypes_adapter"]) if terminal is None else None
    except UnsupportedTandemGenotypesVersion: terminal="UNSUPPORTED_VERSION"; adapter=None
    runs=[]; all_outputs=[]; all_records=[]
    if terminal:
        atomic_write_json(out/"stage-outputs.json",{"record_schema_version":"1.0","outputs":[]}); atomic_write_json(out/"stage-normalized.json",{"record_schema_version":"1.0","evidence_state":"NOT_COMPUTED","records":[],"normalization_warnings":[terminal]}); atomic_write_json(out/"stage-summary.json",{"record_schema_version":"1.0","stage_id":"07_run_tandem_genotypes","status":terminal,"runs":[],"warnings":[terminal],"failure":None})
        if required and not dry_run: raise RuntimeError(f"required tandem-genotypes unavailable: {terminal}")
        return []
    for prep in preparations:
        if prep["status"]!="SUCCEEDED": continue
        seqid=prep["sequence_id"]; alignment=Path(prep["alignment_path"])
        if not alignment.is_file() or __import__('hashlib').sha256(alignment.read_bytes()).hexdigest()!=prep["alignment_file_sha256"]: raise ValueError(f"alignment identity changed: {seqid}")
        final=out/seqid
        if final.exists() and not overwrite: raise FileExistsError(f"caller output exists: {final}")
        work=out/f".{seqid}-work-{uuid4()}"; native=work/"native"; native.mkdir(parents=True); (work/"execution-records").mkdir()
        primary=native/f"{seqid}.tsv"; plan=adapter.plan(tool.resolved_executable,alignment=alignment,repeat_definition=repeat,output=primary,additional_arguments=resources["tandem_genotypes_additional_arguments"])
        recpath=work/"execution-records"/"tandem-genotypes.json"; execute(CommandSpec(f"tandem-genotypes-{seqid}","07_run_tandem_genotypes","TANDEM_GENOTYPES",plan.argv,str(work.resolve()),declared_inputs=(InputDeclaration("alignment",str(alignment)),InputDeclaration("repeat-definition",str(repeat))),declared_outputs=(OutputDeclaration("primary-native",str(primary)),),overwrite=True),tool,recpath,work/"logs")
        produced=sorted(p for p in native.iterdir() if p.is_file()); outputs=[NativeCallerOutput.from_path(p,file_id=f"tandem-genotypes-{seqid}-{i+1}",caller="TANDEM_GENOTYPES",caller_version=tool.detected_version,analysis_source="ASSEMBLED_CONTIG",producer_command_id=f"tandem-genotypes-{seqid}") for i,p in enumerate(produced)]
        status="SUCCEEDED"; warnings=["provisional tandem-genotypes adapter"]
        try: normalized=[normalize_record(r,config=config,caller_version=tool.detected_version,source=outputs[0],associated_sequence_id=seqid) for r in parse_native_tsv(primary)]
        except UnsupportedTandemGenotypesFormat as exc: normalized=[]; status="UNSUPPORTED_FORMAT"; warnings.append(str(exc))
        run={"record_schema_version":"1.0","run_id":f"tandem-genotypes-{seqid}","stage_id":"07_run_tandem_genotypes","caller":"TANDEM_GENOTYPES","caller_version":tool.detected_version,"raw_version_output":tool.raw_version_output,"analysis_source":"ASSEMBLED_CONTIG","associated_sequence_id":seqid,"assignment_state":"DIRECT_SEQUENCE_ASSOCIATION","status":status,"alignment_identity":asdict(file_identity(alignment)),"repeat_definition_identity":asdict(file_identity(repeat)),"command_record_paths":[str(final/"execution-records"/"tandem-genotypes.json")],"native_output_file_ids":[x.file_id for x in outputs],"started_utc":utc_now(),"completed_utc":utc_now(),"warnings":warnings,"failure":None}
        published_outputs=[{**x.to_dict(),"path":str(final/"native"/Path(x.path).name),"associated_sequence_id":seqid} for x in outputs]
        atomic_write_json(work/"run.json",run); atomic_write_json(work/"outputs.json",{"record_schema_version":"1.0","outputs":published_outputs}); atomic_write_json(work/"normalized.json",{"record_schema_version":"1.0","evidence_state":"AVAILABLE" if normalized else "UNSUPPORTED_FORMAT","records":normalized,"normalization_warnings":warnings})
        old=out/f".{seqid}-old-{uuid4()}";
        if final.exists(): os.replace(final,old)
        os.replace(work,final)
        if old.exists(): shutil.rmtree(old)
        # rewrite paths after atomic directory publication
        all_outputs.extend(published_outputs)
        runs.append(run); all_records.extend(normalized)
    status="UNSUPPORTED_FORMAT" if any(r["status"]=="UNSUPPORTED_FORMAT" for r in runs) else "SUCCEEDED"
    atomic_write_json(out/"stage-outputs.json",{"record_schema_version":"1.0","outputs":all_outputs}); atomic_write_json(out/"stage-normalized.json",{"record_schema_version":"1.0","evidence_state":"AVAILABLE" if all_records else "UNSUPPORTED_FORMAT","records":all_records,"normalization_warnings":[]}); atomic_write_json(out/"stage-summary.json",{"record_schema_version":"1.0","stage_id":"07_run_tandem_genotypes","status":status,"runs":runs,"warnings":[],"failure":None}); return runs
