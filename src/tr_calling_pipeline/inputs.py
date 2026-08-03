"""Immutable source-input validation and stage-00 report generation."""
from __future__ import annotations
from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Any

from .config import load_locus_config
from .fasta import read_fasta, validation_report
from .provenance import atomic_write_json, sha256_file
from .sequence_metadata import select_patient_sequences


REFERENCE_SCOPES = {"WHOLE_GENOME_REFERENCE", "LOCAL_LOCUS_REFERENCE", "UNKNOWN_REFERENCE_SCOPE"}


@dataclass(frozen=True)
class ReferenceSequenceSummary:
    record_schema_version: str
    record_id: str
    sequence_length: int


def parse_fai(path: str | Path) -> tuple[ReferenceSequenceSummary, ...]:
    source = Path(path)
    if not source.is_file() or source.stat().st_size == 0:
        raise ValueError(f"reference FASTA index is missing or empty: {source}")
    entries: list[ReferenceSequenceSummary] = []
    for number, line in enumerate(source.read_text(encoding="utf-8").splitlines(), 1):
        fields = line.split("\t")
        if len(fields) < 5 or not fields[0]:
            raise ValueError(f"invalid FAI entry at line {number}")
        try:
            length, offset, line_bases, line_width = map(int, fields[1:5])
        except ValueError as exc:
            raise ValueError(f"invalid numeric FAI field at line {number}") from exc
        if length <= 0 or offset < 0 or line_bases <= 0 or line_width < line_bases:
            raise ValueError(f"invalid FAI dimensions at line {number}")
        entries.append(ReferenceSequenceSummary("1.0", fields[0], length))
    names = [entry.record_id for entry in entries]
    if len(names) != len(set(names)):
        raise ValueError("duplicate reference FASTA index record identifier")
    return tuple(entries)


def validate_reference(fasta: str | Path, fai: str | Path, *, scope: str = "UNKNOWN_REFERENCE_SCOPE", required_contig: str | None = None) -> dict[str, Any]:
    if scope not in REFERENCE_SCOPES:
        raise ValueError(f"unsupported reference_scope: {scope}")
    fasta_records = read_fasta(fasta)
    index = parse_fai(fai)
    fasta_lengths = {r.identifier: len(r.sequence) for r in fasta_records}
    index_lengths = {r.record_id: r.sequence_length for r in index}
    if fasta_lengths != index_lengths:
        raise ValueError("reference FASTA records/lengths do not agree with the configured FAI")
    if required_contig and required_contig not in index_lengths:
        raise ValueError(f"configured locus contig is absent from reference index: {required_contig}")
    return {"schema_version": "1.0", "valid": True, "reference_scope": scope,
        "reference_fasta_path": str(Path(fasta).resolve()), "reference_fasta_sha256": sha256_file(fasta),
        "reference_fasta_index_path": str(Path(fai).resolve()), "reference_fasta_index_sha256": sha256_file(fai),
        "records": [asdict(item) for item in index], "warnings": []}


def _locus_contig(locus: dict[str, Any]) -> str | None:
    locus_data = locus.get("locus", {})
    region = locus_data.get("target_region", {})
    return region.get("contig") or region.get("chromosome") or locus_data.get("chromosome")


def validate_inputs(config: dict[str, Any], output_directory: str | Path) -> dict[str, Any]:
    """Perform Python-only validation; BAM structure belongs to stage 01."""
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    source_paths = [Path(config["inputs"][name]) for name in ("assembly_fasta", "mini_bam", "mini_bam_index", "reference_fasta", "reference_fasta_index")]
    source_paths.append(Path(config["locus_config"]))
    for path in source_paths:
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"required input is missing or empty: {path}")
    assembly = validation_report(config["inputs"]["assembly_fasta"])
    selected, metadata = select_patient_sequences(config["inputs"]["assembly_fasta"], config["inputs"]["assembly_records"])
    locus = load_locus_config(config["locus_config"])
    reference = validate_reference(config["inputs"]["reference_fasta"], config["inputs"]["reference_fasta_index"],
        scope=config["inputs"].get("reference_scope", "UNKNOWN_REFERENCE_SCOPE"), required_contig=_locus_contig(locus))
    metadata_document = {"schema_version": "1.0", "sequences": [item.to_dict() for item in metadata]}
    validated = {"schema_version": "1.0", "valid": True, "division_of_responsibility":
        "Stage 00 validates files, FASTA/config selection, and FAI agreement; stage 01 uses samtools for BAM structure/index validation.",
        "input_files": [{"path": str(path.resolve()), "sha256": sha256_file(path)} for path in source_paths],
        "selected_sequence_count": len(selected), "warnings": []}
    atomic_write_json(output / "assembly-fasta.validation.json", assembly)
    atomic_write_json(output / "reference-fasta.validation.json", reference)
    atomic_write_json(output / "patient-sequences.metadata.json", metadata_document)
    atomic_write_json(output / "validated-inputs.json", validated)
    return validated
