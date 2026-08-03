"""Development-gated LAST preparation primitives.

No LAST release is claimed as laboratory verified.  Command spelling is kept in
one provisional adapter and is never selected without explicit opt-in.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from enum import Enum
from pathlib import Path
import re
import json, os, shutil
from uuid import uuid4

from .provenance import sha256_file
from .tools import Tool, ToolId, detect_version, resolve_tool
from .execution import CommandSpec, InputDeclaration, OutputDeclaration, execute
from .provenance import atomic_write_json, file_identity, utc_now
from .config import load_locus_config


class LastVersionClassification(str, Enum):
    VERIFIED_SUPPORTED = "VERIFIED_SUPPORTED"
    PROVISIONAL_DEVELOPMENT = "PROVISIONAL_DEVELOPMENT"
    UNSUPPORTED = "UNSUPPORTED"
    UNDETERMINED = "UNDETERMINED"


class UnsupportedLastVersion(ValueError): pass
class UnsupportedLastAlignment(ValueError): pass


@dataclass(frozen=True)
class LastCapabilities:
    alignment_format: str = "MAF"
    database_signature_suffixes: tuple[str, ...] = (".prj",)
    lastal_stdout_is_alignment: bool = True


@dataclass(frozen=True)
class LastCommandPlan:
    command_id: str
    argv: tuple[str, ...]
    stdout_path: str | None = None


@dataclass(frozen=True)
class LastAdapter:
    name: str = "last-development-provisional"
    capabilities: LastCapabilities = LastCapabilities()

    def lastdb_plan(self, executable: str, database_prefix: Path, fasta: Path,
                    additional_arguments=()) -> LastCommandPlan:
        return LastCommandPlan("lastdb", (executable, *tuple(additional_arguments), str(database_prefix), str(fasta)))

    def lastal_plan(self, executable: str, database_prefix: Path, reference: Path,
                    alignment: Path, additional_arguments=()) -> LastCommandPlan:
        return LastCommandPlan("lastal", (executable, *tuple(additional_arguments), str(database_prefix), str(reference)), str(alignment))


def classify_version(version: str | None) -> LastVersionClassification:
    if version is None: return LastVersionClassification.UNDETERMINED
    # LAST uses date-like numeric releases. Syntax remains explicitly provisional.
    if re.fullmatch(r"[0-9]{3,5}", version): return LastVersionClassification.PROVISIONAL_DEVELOPMENT
    return LastVersionClassification.UNSUPPORTED


def select_adapter(version: str | None, *, allow_provisional: bool = False) -> LastAdapter:
    kind = classify_version(version)
    if kind is LastVersionClassification.PROVISIONAL_DEVELOPMENT and allow_provisional: return LastAdapter()
    if kind is LastVersionClassification.UNDETERMINED:
        raise UnsupportedLastVersion("LAST version could not be determined; syntax will not be guessed")
    if kind is LastVersionClassification.PROVISIONAL_DEVELOPMENT:
        raise UnsupportedLastVersion("LAST has no laboratory-verified adapter; provisional adapter is disabled")
    raise UnsupportedLastVersion(f"unsupported LAST version: {version}")


def detect_last(executable: str, config_directory: Path, required: bool, tool_id: ToolId) -> Tool:
    if tool_id not in (ToolId.LASTDB, ToolId.LASTAL): raise ValueError("tool_id must be LASTDB or LASTAL")
    tool = resolve_tool(Tool(tool_id, tool_id.value.lower(), executable, required), config_directory)
    return detect_version(tool, pattern=r"(?i)(?:last(?:db|al)?[^0-9]*)?([0-9]{3,5})") if tool.resolved_executable else tool


def write_single_record_fasta(path: str | Path, sequence_id: str, sequence: str, width: int = 80) -> str:
    """Write deterministic FASTA without altering sequence letters."""
    if not sequence or "\n" in sequence or "\r" in sequence: raise ValueError("sequence must be non-empty and newline-free")
    target = Path(path); target.parent.mkdir(parents=True, exist_ok=True)
    text = f">{sequence_id}\n" + "\n".join(sequence[i:i+width] for i in range(0, len(sequence), width)) + "\n"
    target.write_text(text, encoding="ascii", newline="\n")
    return sha256_file(target)


def validate_maf(path: str | Path, expected_query_id: str) -> tuple[str, ...]:
    source = Path(path)
    if not source.is_file(): raise UnsupportedLastAlignment("alignment is missing")
    if source.stat().st_size == 0: raise UnsupportedLastAlignment("alignment is empty")
    lines = source.read_text(encoding="utf-8").splitlines()
    if not lines or not lines[0].startswith("##maf"): raise UnsupportedLastAlignment("unsupported alignment format")
    blocks = [i for i, line in enumerate(lines) if line.startswith("a ") or line == "a"]
    if not blocks: raise UnsupportedLastAlignment("MAF has no alignment block")
    identifiers = tuple(line.split()[1] for line in lines if line.startswith("s ") and len(line.split()) >= 7)
    if expected_query_id not in identifiers: raise UnsupportedLastAlignment("expected patient sequence is absent from MAF")
    return identifiers

def _fasta(path: Path) -> dict[str, str]:
    records={}; name=None; chunks=[]
    for line in path.read_text(encoding="ascii").splitlines():
        if line.startswith(">"):
            if name is not None: records[name]="".join(chunks)
            name=line[1:].split()[0]; chunks=[]
        elif name is None: raise ValueError("FASTA sequence precedes header")
        else: chunks.append(line)
    if name is not None: records[name]="".join(chunks)
    return records

def _resource(config):
    locus=Path(config["locus_config"]); value=load_locus_config(locus)["caller_resources"]["tandem_genotypes"]["repeat_definition"]
    if value is None: return None
    p=Path(value).expanduser(); return p if p.is_absolute() else (locus.parent/p).resolve()

def prepare_tandem_genotypes_stage(config, root: Path, config_directory: Path, *, overwrite=False, dry_run=False):
    """Prepare one isolated LAST database and MAF per authoritative sequence."""
    out=root/"07_tandem_genotypes_preparation"; out.mkdir(parents=True,exist_ok=True)
    summary_path=out/"preparation-summary.json"; align_registry=out/"alignment-inputs.json"; metadata_registry=out/"alignment-metadata.json"
    settings=load_locus_config(config["locus_config"])["caller_resources"]["tandem_genotypes"]
    tools=[]
    for tid,key in ((ToolId.LASTDB,"lastdb"),(ToolId.LASTAL,"lastal")):
        tc=config.get("tools",{}).get(key,{}); tools.append(detect_last(tc.get("executable",key),config_directory,tc.get("required",False),tid))
    repeat=_resource(config); status=None; warning=None; dry_warnings=[]
    if repeat is None or not repeat.is_file():
        status="INPUT_MISSING"; warning="explicit repeat definition is missing"
        dry_warnings.append(warning)
    for tool in tools:
        if not tool.resolved_executable:
            dry_warnings.append(f"configured {tool.tool_id.value} executable is unavailable")
    if status is None and any(not t.resolved_executable for t in tools): status="TOOL_MISSING"; warning="configured LAST executable is unavailable"
    adapters=[]
    if status is None:
        try: adapters=[select_adapter(t.detected_version,allow_provisional=settings["allow_provisional_last_adapter"]) for t in tools]
        except UnsupportedLastVersion as exc: status="UNSUPPORTED_VERSION"; warning=str(exc)
    if dry_run:
        # Probe and report blockers, but a dry run is always truthfully a plan.
        if all(t.resolved_executable for t in tools):
            for tool in tools:
                try: select_adapter(tool.detected_version,allow_provisional=settings["allow_provisional_last_adapter"])
                except UnsupportedLastVersion as exc: dry_warnings.append(f"{tool.tool_id.value}: {exc}")
        status="DRY_RUN"; warning="; ".join(dict.fromkeys(dry_warnings)) or "Dry run: LAST commands were planned but not executed"
    records=[]
    if status:
        patient_metadata=root/"03_assembly_alignment"/"patient-sequences.metadata.json"
        metadata_records=json.loads(patient_metadata.read_text()).get("sequences",[]) if patient_metadata.is_file() else []
        if dry_run and not patient_metadata.is_file():
            warning += "; authoritative patient metadata is unavailable"
        for m in metadata_records:
            records.append({"record_schema_version":"1.0","preparation_id":f"last-{m['sequence_id']}","stage_id":"06_prepare_tandem_genotypes","sequence_id":m["sequence_id"],"source_fasta_record_id":m["source_fasta_record_id"],"sequence_sha256":m["sequence_sha256"],"status":status,"lastdb_version":tools[0].detected_version,"lastal_version":tools[1].detected_version,"input_file_ids":[],"database_file_ids":[],"alignment_file_id":None,"alignment_file_sha256":None,"alignment_path":"","alignment_format":None,"coordinate_space":"UNKNOWN_COORDINATE_SPACE","repeat_definition_identity":asdict(file_identity(repeat)) if repeat and repeat.is_file() else None,"command_record_paths":[],"started_utc":utc_now(),"completed_utc":utc_now(),"warnings":[warning] if warning else [],"failure":None})
        warnings=[warning] if warning else []
        atomic_write_json(align_registry,{"record_schema_version":"1.0","alignments":[]}); atomic_write_json(metadata_registry,{"record_schema_version":"1.0","records":records}); atomic_write_json(summary_path,{"record_schema_version":"1.0","stage_id":"06_prepare_tandem_genotypes","status":status,"preparations":records,"warnings":warnings,"failure":None})
        if not dry_run and any(t.required for t in tools): raise RuntimeError(warning)
        return records
    fasta=root/"03_assembly_alignment"/"patient-sequences.fasta"; metadata_path=root/"03_assembly_alignment"/"patient-sequences.metadata.json"
    sequences=_fasta(fasta); metadata=json.loads(metadata_path.read_text())["sequences"]
    for m in metadata:
        seqid=m["sequence_id"]; seq=sequences.get(seqid)
        if seq is None or __import__('hashlib').sha256(seq.encode('ascii')).hexdigest()!=m["sequence_sha256"]: raise ValueError(f"patient sequence checksum mismatch: {seqid}")
        final=out/seqid
        if final.exists() and not overwrite: raise FileExistsError(f"preparation exists: {final}")
        work=out/f".{seqid}-work-{uuid4()}"; (work/"input").mkdir(parents=True); (work/"lastdb").mkdir(); (work/"alignment").mkdir(); (work/"execution-records").mkdir()
        input_fasta=work/"input"/f"{seqid}.fasta"; input_sha=write_single_record_fasta(input_fasta,seqid,seq)
        prefix=work/"lastdb"/"database"; dbplan=adapters[0].lastdb_plan(tools[0].resolved_executable,prefix,input_fasta,settings["lastdb_additional_arguments"])
        dbrec=work/"execution-records"/"lastdb.json"; execute(CommandSpec(f"lastdb-{seqid}","06_prepare_tandem_genotypes","LASTDB",dbplan.argv,str(work.resolve()),declared_inputs=(InputDeclaration("patient-fasta",str(input_fasta)),),overwrite=True),tools[0],dbrec,work/"logs")
        dbfiles=sorted(p for p in (work/"lastdb").iterdir() if p.is_file())
        if not dbfiles or not any(p.suffix in adapters[0].capabilities.database_signature_suffixes for p in dbfiles): raise RuntimeError("LASTDB signature output missing")
        alignment=work/"alignment"/f"{seqid}.maf"; alplan=adapters[1].lastal_plan(tools[1].resolved_executable,prefix,Path(config["inputs"]["reference_fasta"]),alignment,settings["lastal_additional_arguments"])
        alrec=work/"execution-records"/"lastal.json"; execute(CommandSpec(f"lastal-{seqid}","06_prepare_tandem_genotypes","LASTAL",alplan.argv,str(work.resolve()),declared_inputs=(InputDeclaration("database-signature",str(dbfiles[0])),InputDeclaration("reference",config["inputs"]["reference_fasta"])),overwrite=True),tools[1],alrec,work/"logs")
        stdout=work/"logs"/"06_prepare_tandem_genotypes"/f"lastal-{seqid}.stdout.log"; shutil.copyfile(stdout,alignment); validate_maf(alignment,seqid)
        old=out/f".{seqid}-old-{uuid4()}";
        if final.exists(): os.replace(final,old)
        os.replace(work,final)
        if old.exists(): shutil.rmtree(old)
        dbfiles=sorted((final/"lastdb").iterdir()); alignment=final/"alignment"/f"{seqid}.maf"
        rec={"record_schema_version":"1.0","preparation_id":f"last-{seqid}","stage_id":"06_prepare_tandem_genotypes","sequence_id":seqid,"source_fasta_record_id":m["source_fasta_record_id"],"sequence_sha256":m["sequence_sha256"],"input_fasta_sha256":input_sha,"status":"SUCCEEDED","lastdb_version":tools[0].detected_version,"lastal_version":tools[1].detected_version,"input_file_ids":[asdict(file_identity(fasta)),asdict(file_identity(metadata_path))],"database_file_ids":[asdict(file_identity(p)) for p in dbfiles],"alignment_file_id":f"last-alignment-{seqid}","alignment_file_sha256":sha256_file(alignment),"alignment_path":str(alignment),"alignment_format":"MAF","coordinate_space":"ALIGNMENT_COORDINATES","repeat_definition_identity":asdict(file_identity(repeat)),"command_record_paths":[str(final/"execution-records"/"lastdb.json"),str(final/"execution-records"/"lastal.json")],"started_utc":utc_now(),"completed_utc":utc_now(),"warnings":["provisional LAST adapter"],"failure":None}; atomic_write_json(final/"alignment.metadata.json",rec); records.append(rec)
    atomic_write_json(align_registry,{"record_schema_version":"1.0","alignments":[{"sequence_id":r["sequence_id"],"path":r["alignment_path"],"sha256":r["alignment_file_sha256"]} for r in records]}); atomic_write_json(metadata_registry,{"record_schema_version":"1.0","records":records}); atomic_write_json(summary_path,{"record_schema_version":"1.0","stage_id":"06_prepare_tandem_genotypes","status":"SUCCEEDED","preparations":records,"warnings":[],"failure":None}); return records
