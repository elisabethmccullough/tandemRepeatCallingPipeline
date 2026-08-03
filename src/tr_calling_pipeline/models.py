"""Versioned, immutable shared-contract models and stable vocabularies."""

from __future__ import annotations

from dataclasses import dataclass, field
from decimal import Decimal
from enum import StrEnum
from types import MappingProxyType
from typing import Any, Mapping


class SequenceRole(StrEnum):
    PATIENT_HAPLOTYPE = "PATIENT_HAPLOTYPE"
    PATIENT_CONSENSUS = "PATIENT_CONSENSUS"
    REFERENCE_SEQUENCE = "REFERENCE_SEQUENCE"
    COMPARISON_SEQUENCE = "COMPARISON_SEQUENCE"


class EvidenceState(StrEnum):
    AVAILABLE = "AVAILABLE"
    NOT_APPLICABLE = "NOT_APPLICABLE"
    INPUT_MISSING = "INPUT_MISSING"
    COMPUTATION_FAILED = "COMPUTATION_FAILED"
    AMBIGUOUS = "AMBIGUOUS"
    NOT_COMPUTED = "NOT_COMPUTED"
    UNSUPPORTED_FORMAT = "UNSUPPORTED_FORMAT"


class AnalysisSource(StrEnum):
    RAW_READS = "RAW_READS"
    ASSEMBLED_CONTIG = "ASSEMBLED_CONTIG"


class AssignmentState(StrEnum):
    DIRECT_SEQUENCE_ASSOCIATION = "DIRECT_SEQUENCE_ASSOCIATION"
    UNASSIGNED = "UNASSIGNED"
    ASSIGNMENT_NOT_APPLICABLE = "ASSIGNMENT_NOT_APPLICABLE"


@dataclass(frozen=True, slots=True)
class AssemblyRecord:
    record_id: str
    sequence_id: str
    display_label: str
    sequence_role: SequenceRole


@dataclass(frozen=True, slots=True)
class RunConfiguration:
    """Typed identity and assembly portion of a validated run configuration."""
    schema_version: str
    case_id: str
    subject_id: str
    sample_id: str
    locus_id: str
    assembly_records: tuple[AssemblyRecord, ...]


@dataclass(frozen=True, slots=True)
class LocusConfigurationIdentity:
    schema_version: str
    locus_config_id: str
    locus_config_version: str
    release_status: str
    locus_id: str


@dataclass(frozen=True, slots=True)
class NativeCallerOutput:
    caller: str
    caller_version: str
    analysis_source: AnalysisSource
    file_id: str
    status: EvidenceState


@dataclass(frozen=True, slots=True)
class FileInventoryItem:
    file_id: str
    role: str
    path: str
    required: bool
    status: EvidenceState
    sha256: str | None
    size_bytes: int | None
    media_type: str
    producer_stage: str


@dataclass(frozen=True, slots=True)
class CaseManifestIdentity:
    schema_version: str
    case_id: str
    subject_id: str
    package_created_utc: str


def _freeze(value: Any) -> Any:
    if isinstance(value, dict):
        return MappingProxyType({key: _freeze(item) for key, item in value.items()})
    if isinstance(value, list):
        return tuple(_freeze(item) for item in value)
    return value


@dataclass(frozen=True, slots=True)
class CallerEvidence:
    """One normalized evidence record; strings and Decimal values are lossless."""

    record_schema_version: str
    record_id: str
    case_id: str
    subject_id: str
    sample_id: str
    locus_id: str
    caller: str
    caller_version: str
    analysis_source: AnalysisSource
    source_file_id: str
    source_file_sha256: str
    native_record_identifier: str | None
    native_allele_identifier: str | None
    associated_sequence_id: str | None
    assignment_state: AssignmentState
    reference_build: str
    chromosome: str
    start: int | None
    end: int | None
    coordinate_convention: str
    reported_motif: str | None
    reported_motif_chain: tuple[Any, ...] | None
    reported_repeat_count: str | Decimal | None
    reported_repeat_length_bp: str | Decimal | None
    supporting_reads: int | None
    total_spanning_reads: int | None
    quality_state: EvidenceState
    raw_fields: Mapping[str, Any] = field(default_factory=dict)
    normalization_warnings: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "raw_fields", _freeze(dict(self.raw_fields)))
        if self.reported_motif_chain is not None:
            object.__setattr__(self, "reported_motif_chain", tuple(_freeze(x) for x in self.reported_motif_chain))
        object.__setattr__(self, "normalization_warnings", tuple(self.normalization_warnings))
        if self.analysis_source is AnalysisSource.RAW_READS and self.assignment_state is AssignmentState.DIRECT_SEQUENCE_ASSOCIATION:
            raise ValueError("inferred raw-read sequence assignment is not supported")
