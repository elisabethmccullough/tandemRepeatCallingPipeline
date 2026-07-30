"""Common result schema and output writers; caller parsers will be added later."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Iterable, Mapping, Any

NORMALIZED_COLUMNS = [
    "sample_id", "locus_id", "caller", "analysis_source", "reference_build",
    "chrom", "start", "end", "reference_motif", "reported_motifs", "allele_id",
    "estimated_repeat_count", "estimated_repeat_length_bp", "repeat_structure",
    "supporting_reads", "total_spanning_reads", "strand_forward_reads",
    "strand_reverse_reads", "haplotype", "quality_status", "native_output_path",
    "normalization_notes",
]
ALLOWED_ANALYSIS_SOURCES = frozenset({"raw_reads", "assembled_contig"})
WARNING_COLUMNS = ["caller", "native_output_path", "warning"]


def validate_record(record: Mapping[str, Any]) -> None:
    source = record.get("analysis_source")
    if source not in ALLOWED_ANALYSIS_SOURCES:
        raise ValueError(f"Unsupported analysis_source: {source!r}")


def write_normalized_outputs(
    output_dir: str | Path,
    records: Iterable[Mapping[str, Any]] = (),
    warnings: Iterable[Mapping[str, Any]] = (),
) -> dict[str, Path]:
    """Write deterministic TSV/JSON outputs, including headers when empty."""
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    rows = list(records)
    for row in rows:
        validate_record(row)
    paths = {
        "tsv": destination / "caller_results.tsv",
        "json": destination / "caller_results.json",
        "warnings": destination / "normalization_warnings.tsv",
    }
    with paths["tsv"].open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=NORMALIZED_COLUMNS, extrasaction="ignore", delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    paths["json"].write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
    with paths["warnings"].open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=WARNING_COLUMNS, extrasaction="ignore", delimiter="\t")
        writer.writeheader()
        writer.writerows(warnings)
    return paths
