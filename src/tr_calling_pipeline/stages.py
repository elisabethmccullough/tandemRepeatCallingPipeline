"""Canonical stage registry and conservative resume evaluation."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Mapping, Sequence


class StageStatus(str, Enum):
    NOT_STARTED = "NOT_STARTED"
    PLANNED = "PLANNED"
    RUNNING = "RUNNING"
    SUCCEEDED = "SUCCEEDED"
    FAILED = "FAILED"
    SKIPPED = "SKIPPED"
    DRY_RUN = "DRY_RUN"
    INVALIDATED = "INVALIDATED"


class ResumeReason(str, Enum):
    RESUME_ALLOWED = "RESUME_ALLOWED"
    RESUME_DISABLED = "RESUME_DISABLED"
    NO_PRIOR_RECORD = "NO_PRIOR_RECORD"
    PRIOR_STAGE_FAILED = "PRIOR_STAGE_FAILED"
    INPUT_CHANGED = "INPUT_CHANGED"
    OUTPUT_MISSING = "OUTPUT_MISSING"
    OUTPUT_CHANGED = "OUTPUT_CHANGED"
    CONFIGURATION_CHANGED = "CONFIGURATION_CHANGED"
    TOOL_CHANGED = "TOOL_CHANGED"
    TOOL_VERSION_CHANGED = "TOOL_VERSION_CHANGED"
    EXECUTION_MODE_CHANGED = "EXECUTION_MODE_CHANGED"
    CONTAINER_CHANGED = "CONTAINER_CHANGED"
    UNSUPPORTED_RECORD_VERSION = "UNSUPPORTED_RECORD_VERSION"


@dataclass(frozen=True)
class StageDefinition:
    stage_id: str
    order: int
    display_name: str
    description: str
    required_input_roles: tuple[str, ...] = ()
    expected_output_roles: tuple[str, ...] = ()
    required_tools: tuple[str, ...] = ()
    optional_tools: tuple[str, ...] = ()
    supports_dry_run: bool = True
    skippable: bool = True
    resume_supported: bool = True


def _stage(stage_id: str, order: int, inputs=(), outputs=(), required=(), optional=()) -> StageDefinition:
    return StageDefinition(
        stage_id=stage_id,
        order=order,
        display_name=stage_id.replace("_", " ").title(),
        description=f"Scaffold workflow stage {stage_id}.",
        required_input_roles=tuple(inputs),
        expected_output_roles=tuple(outputs),
        required_tools=tuple(required),
        optional_tools=tuple(optional),
    )


STAGES = (
    _stage("00_validate_inputs", 0, ("ASSEMBLY_FASTA", "MINI_BAM", "MINI_BAM_INDEX", "REFERENCE_FASTA", "REFERENCE_FASTA_INDEX", "LOCUS_CONFIG"), ("VALIDATED_INPUTS", "PATIENT_SEQUENCE_METADATA"), ("PYTHON",)),
    _stage("01_prepare_bam", 1, ("MINI_BAM", "MINI_BAM_INDEX"), ("PREPARED_MINI_BAM", "PREPARED_MINI_BAM_INDEX", "BAM_VALIDATION_REPORT", "BAM_FLAGSTAT", "BAM_IDXSTATS"), ("SAMTOOLS",)),
    _stage("02_align_assembly", 2, ("ASSEMBLY_FASTA", "REFERENCE_FASTA", "REFERENCE_FASTA_INDEX", "PATIENT_SEQUENCE_METADATA"), ("ALIGNED_ASSEMBLY", "ALIGNED_ASSEMBLY_INDEX", "ASSEMBLY_ALIGNMENT_SUMMARY", "ASSEMBLY_RECORD_MAPPINGS", "PACKAGE_PATIENT_FASTA", "PACKAGE_PATIENT_METADATA"), ("MINIMAP2", "SAMTOOLS")),
    _stage("03_run_vamos_read", 3, ("PREPARED_MINI_BAM", "PREPARED_MINI_BAM_INDEX", "REFERENCE_FASTA", "REFERENCE_FASTA_INDEX", "LOCUS_CONFIG"), ("VAMOS_READ_NATIVE_OUTPUTS", "VAMOS_READ_RUN_METADATA", "VAMOS_READ_NORMALIZED_EVIDENCE"), optional=("VAMOS",)),
    _stage("04_run_vamos_contig", 4, ("PACKAGE_PATIENT_FASTA", "PACKAGE_PATIENT_METADATA", "ASSEMBLY_RECORD_MAPPINGS", "REFERENCE_FASTA", "REFERENCE_FASTA_INDEX", "LOCUS_CONFIG"), ("VAMOS_CONTIG_NATIVE_OUTPUTS", "VAMOS_CONTIG_RUN_METADATA", "VAMOS_CONTIG_NORMALIZED_EVIDENCE"), optional=("VAMOS",)),
    _stage("05_run_straglr", 5, ("PREPARED_MINI_BAM", "PREPARED_MINI_BAM_INDEX", "REFERENCE_FASTA", "REFERENCE_FASTA_INDEX", "LOCUS_CONFIG"), ("STRAGLR_NATIVE_OUTPUTS", "STRAGLR_RUN_METADATA", "STRAGLR_NORMALIZED_EVIDENCE"), optional=("STRAGLR",)),
    _stage("06_prepare_tandem_genotypes", 6, ("PACKAGE_PATIENT_FASTA", "PACKAGE_PATIENT_METADATA", "ASSEMBLY_RECORD_MAPPINGS", "REFERENCE_FASTA", "REFERENCE_FASTA_INDEX", "LOCUS_CONFIG"), ("TANDEM_GENOTYPES_ALIGNMENT_INPUTS", "TANDEM_GENOTYPES_ALIGNMENT_METADATA", "TANDEM_GENOTYPES_PREPARATION_SUMMARY"), optional=("LASTDB", "LASTAL")),
    _stage("07_run_tandem_genotypes", 7, ("TANDEM_GENOTYPES_ALIGNMENT_INPUTS", "TANDEM_GENOTYPES_ALIGNMENT_METADATA", "PACKAGE_PATIENT_METADATA", "LOCUS_CONFIG"), ("TANDEM_GENOTYPES_NATIVE_OUTPUTS", "TANDEM_GENOTYPES_RUN_METADATA", "TANDEM_GENOTYPES_NORMALIZED_EVIDENCE"), optional=("TANDEM_GENOTYPES",)),
    _stage("08_normalize_outputs", 8, (
        "VAMOS_READ_NATIVE_OUTPUTS", "VAMOS_READ_RUN_METADATA", "VAMOS_READ_NORMALIZED_EVIDENCE",
        "VAMOS_CONTIG_NATIVE_OUTPUTS", "VAMOS_CONTIG_RUN_METADATA", "VAMOS_CONTIG_NORMALIZED_EVIDENCE",
        "STRAGLR_NATIVE_OUTPUTS", "STRAGLR_RUN_METADATA", "STRAGLR_NORMALIZED_EVIDENCE",
        "TANDEM_GENOTYPES_NATIVE_OUTPUTS", "TANDEM_GENOTYPES_RUN_METADATA", "TANDEM_GENOTYPES_NORMALIZED_EVIDENCE",
        "PACKAGE_PATIENT_FASTA", "PACKAGE_PATIENT_METADATA", "LOCUS_CONFIG"),
        ("UNIFIED_NORMALIZED_EVIDENCE", "UNIFIED_EVIDENCE_SUMMARY", "UNIFIED_SOURCE_REGISTRY", "UNIFIED_VALIDATION_REPORT"), ("PYTHON",)),
    _stage("09_build_case_package", 9, ("NORMALIZED_EVIDENCE", "ASSEMBLY_FASTA"), ("CASE_PACKAGE",), ("PIPELINE",)),
    _stage("10_validate_case_package", 10, ("CASE_PACKAGE",), ("VALIDATED_CASE_PACKAGE",), ("PIPELINE",)),
)

if len({stage.stage_id for stage in STAGES}) != len(STAGES) or len({stage.order for stage in STAGES}) != len(STAGES):
    raise RuntimeError("stage IDs and orders must be unique")


def select_stages(start: str | None = None, stop: str | None = None) -> tuple[StageDefinition, ...]:
    index = {stage.stage_id: position for position, stage in enumerate(STAGES)}
    start = start or STAGES[0].stage_id
    stop = stop or STAGES[-1].stage_id
    if start not in index:
        raise ValueError(f"unknown start stage: {start}")
    if stop not in index:
        raise ValueError(f"unknown stop stage: {stop}")
    if index[start] > index[stop]:
        raise ValueError("start stage follows stop stage")
    return STAGES[index[start] : index[stop] + 1]


def _by_path(identities: Sequence[Mapping[str, object]]) -> dict[str, Mapping[str, object]]:
    return {str(Path(str(identity["path"])).resolve()): identity for identity in identities}


def resume_eligibility(
    prior: Mapping[str, object] | None,
    configuration_digest: str,
    input_identities: list[Mapping[str, object]],
    output_identities: list[Mapping[str, object]],
    tool_identities: list[Mapping[str, object]],
) -> ResumeReason:
    if not prior:
        return ResumeReason.NO_PRIOR_RECORD
    if prior.get("record_schema_version") != "1.0":
        return ResumeReason.UNSUPPORTED_RECORD_VERSION
    if prior.get("status") != StageStatus.SUCCEEDED.value:
        return ResumeReason.PRIOR_STAGE_FAILED
    if prior.get("configuration_digest") != configuration_digest:
        return ResumeReason.CONFIGURATION_CHANGED

    current_inputs = _by_path(input_identities)
    for expected in prior.get("input_file_identities", []):
        current = current_inputs.get(str(Path(str(expected["path"])).resolve()))
        if current != expected:
            return ResumeReason.INPUT_CHANGED

    current_outputs = _by_path(output_identities)
    for expected in prior.get("output_file_identities", []):
        path = str(Path(str(expected["path"])).resolve())
        if path not in current_outputs:
            return ResumeReason.OUTPUT_MISSING
        if current_outputs[path] != expected:
            return ResumeReason.OUTPUT_CHANGED

    old_tools = prior.get("tool_identities", [])
    identity_keys = ("tool_id", "resolved_executable")
    if [{key: tool.get(key) for key in identity_keys} for tool in old_tools] != [
        {key: tool.get(key) for key in identity_keys} for tool in tool_identities
    ]:
        return ResumeReason.TOOL_CHANGED
    if [tool.get("detected_version") for tool in old_tools] != [tool.get("detected_version") for tool in tool_identities]:
        return ResumeReason.TOOL_VERSION_CHANGED
    if [tool.get("execution_mode") for tool in old_tools] != [tool.get("execution_mode") for tool in tool_identities]:
        return ResumeReason.EXECUTION_MODE_CHANGED
    if [tool.get("container_digest") for tool in old_tools] != [tool.get("container_digest") for tool in tool_identities]:
        return ResumeReason.CONTAINER_CHANGED
    return ResumeReason.RESUME_ALLOWED
