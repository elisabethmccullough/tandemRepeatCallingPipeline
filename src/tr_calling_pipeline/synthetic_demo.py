"""Offline synthetic portability demonstration using production package functions."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import shutil

from .case_package import build_case_package
from .case_package_validation import validate_case_package
from .normalization import run_normalization_stage

def run_demo(output: str | Path) -> dict[str, object]:
    """Build, move, and validate a non-clinical package without external tools."""
    output = Path(output).resolve()
    if output.exists(): raise FileExistsError(f"demo output already exists: {output}")
    source = output / "temporary synthetic run"
    patient = source / "03_assembly_alignment"; patient.mkdir(parents=True)
    records = [("synthetic-sequence-1", "ACGTACGT"), ("synthetic-sequence-2", "TGCATGCA")]
    (patient / "patient-sequences.fasta").write_text("".join(f">{name}\n{seq}\n" for name, seq in records), encoding="utf-8")
    sequences = []
    for index, (name, sequence) in enumerate(records, 1):
        sequences.append({"record_schema_version": "1.0", "sequence_id": name, "sequence_role": "PATIENT_HAPLOTYPE",
            "display_label": f"Synthetic sequence {index}", "source_fasta_record_id": f"SYNTHETIC_{index}",
            "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(), "sequence_length": len(sequence),
            "original_orientation": "AS_STORED", "reference_alignment_strand": "NOT_ALIGNED",
            "display_orientation": "UNRESOLVED", "reverse_complement_required": None, "source_coordinates": None,
            "mapping_status": "NOT_ALIGNED", "warnings": ["SYNTHETIC_NON_CLINICAL"]})
    (patient / "patient-sequences.metadata.json").write_text(json.dumps({"schema_version": "1.0", "sequences": sequences}), encoding="utf-8")
    (patient / "assembly_record_mappings.json").write_text(json.dumps({"schema_version": "1.0", "coordinate_convention": "zero_based_half_open", "mappings": []}), encoding="utf-8")
    control = output / "local synthetic fixtures"; control.mkdir(parents=True)
    locus = control / "locus.yaml"; locus.write_text('schema_version: "1.0"\n', encoding="utf-8")
    run_config = control / "run.yaml"; run_config.write_text('schema_version: "1.0"\n', encoding="utf-8")
    identity = {"case_id": "SYNTHETIC-DEMO", "subject_id": "SYNTHETIC-SUBJECT", "sample_id": "SYNTHETIC-SAMPLE", "locus_id": "SYNTHETIC-LOCUS"}
    run_normalization_stage({"run": identity, "locus_config": str(locus)}, source)
    built = output / "built package"
    config = {"run": identity, "locus_config": str(locus), "_config_path": str(run_config), "case_package": {
        "package_root": str(built), "include_native_outputs": True, "include_prepared_mini_bam": False,
        "include_alignment_artifacts": False, "include_command_records": False}}
    build_case_package(config, source)
    before_move = validate_case_package(built)
    moved = output / "portable synthetic case package"; shutil.move(built, moved)
    after_move = validate_case_package(moved)
    shutil.rmtree(source); shutil.rmtree(control)
    after_removal = validate_case_package(moved, write_report=True)
    return {"synthetic": True, "non_clinical": True, "external_tools_executed": False,
            "tool_notice": "No real or fake caller executable was run; this demo verifies contract/package portability only.",
            "package": str(moved), "caller_records_by_source": {}, "valid_before_move": before_move["valid"],
            "valid_after_move": after_move["valid"], "valid_after_source_removal": after_removal["valid"]}
