"""Stable GUI-facing patient sequence contracts."""

from __future__ import annotations
from dataclasses import asdict, dataclass, replace
from enum import Enum
from pathlib import Path
from typing import Any, Iterable

from .fasta import FastaRecord, FastaValidationError, read_fasta, sequence_sha256
from .provenance import sha256_file


class OriginalOrientation(str, Enum):
    AS_STORED = "AS_STORED"


class AlignmentStrand(str, Enum):
    FORWARD = "FORWARD"
    REVERSE = "REVERSE"
    UNKNOWN = "UNKNOWN"
    NOT_ALIGNED = "NOT_ALIGNED"


class DisplayOrientation(str, Enum):
    REFERENCE_FORWARD = "REFERENCE_FORWARD"
    UNRESOLVED = "UNRESOLVED"


class MappingStatus(str, Enum):
    UNIQUE_PRIMARY = "UNIQUE_PRIMARY"
    MULTIPLE_PRIMARY = "MULTIPLE_PRIMARY"
    UNMAPPED = "UNMAPPED"
    AMBIGUOUS = "AMBIGUOUS"
    TRUNCATED = "TRUNCATED"
    SECONDARY_ONLY = "SECONDARY_ONLY"
    SUPPLEMENTARY_ONLY = "SUPPLEMENTARY_ONLY"
    NOT_ALIGNED = "NOT_ALIGNED"


@dataclass(frozen=True)
class PatientSequenceMetadata:
    record_schema_version: str
    sequence_id: str
    sequence_role: str
    display_label: str
    source_fasta_path: str
    source_fasta_file_sha256: str
    source_fasta_record_id: str
    sequence_sha256: str
    sequence_length: int
    original_orientation: str = OriginalOrientation.AS_STORED.value
    reference_alignment_strand: str = AlignmentStrand.UNKNOWN.value
    display_orientation: str = DisplayOrientation.UNRESOLVED.value
    reverse_complement_required: bool | None = None
    source_coordinates: dict[str, Any] | None = None
    mapping_status: str = MappingStatus.NOT_ALIGNED.value
    warnings: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        value = asdict(self)
        value["warnings"] = list(self.warnings)
        return value


def select_patient_sequences(path: str | Path, configured: Iterable[dict[str, Any]]) -> tuple[tuple[FastaRecord, ...], tuple[PatientSequenceMetadata, ...]]:
    source = Path(path).resolve()
    by_id = {record.identifier: record for record in read_fasta(source)}
    definitions = list(configured)
    record_ids = [item["record_id"] for item in definitions]
    sequence_ids = [item["sequence_id"] for item in definitions]
    if len(record_ids) != len(set(record_ids)):
        raise FastaValidationError("configured assembly record_id values must be unique")
    if len(sequence_ids) != len(set(sequence_ids)):
        raise FastaValidationError("configured sequence_id values must be unique")
    allowed = {"PATIENT_HAPLOTYPE", "PATIENT_CONSENSUS"}
    missing = [item for item in record_ids if item not in by_id]
    if missing:
        raise FastaValidationError(f"configured assembly record(s) not found: {', '.join(missing)}")
    digest = sha256_file(source)
    selected: list[FastaRecord] = []
    metadata: list[PatientSequenceMetadata] = []
    for item in definitions:
        if item["sequence_role"] not in allowed:
            raise FastaValidationError(f"unsupported patient sequence role: {item['sequence_role']}")
        record = by_id[item["record_id"]]
        selected.append(FastaRecord(item["sequence_id"], "", record.sequence))
        metadata.append(PatientSequenceMetadata(
            "1.0", item["sequence_id"], item["sequence_role"], item["display_label"], str(source),
            digest, record.identifier, sequence_sha256(record.sequence), len(record.sequence),
        ))
    return tuple(selected), tuple(metadata)


def with_mapping(metadata: PatientSequenceMetadata, *, strand: str, status: str, coordinates: dict[str, Any] | None, warnings: tuple[str, ...] = ()) -> PatientSequenceMetadata:
    resolved = status == MappingStatus.UNIQUE_PRIMARY.value and strand in {"FORWARD", "REVERSE"}
    return replace(metadata, reference_alignment_strand=strand,
        display_orientation=DisplayOrientation.REFERENCE_FORWARD.value if resolved else DisplayOrientation.UNRESOLVED.value,
        reverse_complement_required=(strand == AlignmentStrand.REVERSE.value) if resolved else None,
        source_coordinates=coordinates, mapping_status=status, warnings=warnings)
