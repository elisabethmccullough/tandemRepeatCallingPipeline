"""Lossless parser/normalizer for the provisional tandem-genotypes TSV contract."""
from __future__ import annotations
import csv, re
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any
from ..caller_outputs import NativeCallerOutput

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
