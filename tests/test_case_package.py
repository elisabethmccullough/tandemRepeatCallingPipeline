import hashlib
import json
import shutil
from pathlib import Path

import pytest

from tr_calling_pipeline.case_package import _safe_relative, build_case_package, reject_secret_bearing_config
from tr_calling_pipeline.case_package_validation import validate_case_package
from tr_calling_pipeline.config import ConfigurationError
from tr_calling_pipeline.normalization import run_normalization_stage
from tr_calling_pipeline.provenance import sha256_file


def test_secret_like_configuration_is_rejected(tmp_path):
    config = tmp_path / "unsafe.yaml"
    config.write_text("service:\n  api_key: synthetic-do-not-package\n", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="secret-like"):
        reject_secret_bearing_config(config)


@pytest.mark.parametrize("value", ["../outside", "/absolute", "C:\\windows\\file"])
def test_package_relative_path_security(value):
    with pytest.raises(ConfigurationError):
        _safe_relative(value)


def test_independent_validator_reports_missing_manifest(tmp_path):
    report = validate_case_package(tmp_path)
    assert report["valid"] is False
    assert report["error_count"] > 0
    assert any(issue["code"] == "SCHEMA_VALIDATION_ERROR" for issue in report["issues"])


def test_unexpected_file_is_an_error(tmp_path):
    # A deliberately incomplete package still demonstrates the strict extra-file policy.
    (tmp_path / "extra.txt").write_text("unexpected", encoding="utf-8")
    report = validate_case_package(tmp_path)
    assert any(issue["code"] == "UNEXPECTED_UNREGISTERED_FILE" for issue in report["issues"])


def _built_package(tmp_path):
    run = tmp_path / "source run"
    patient = run / "03_assembly_alignment"
    patient.mkdir(parents=True)
    sequence = "ACGT"
    (patient / "patient-sequences.fasta").write_text(f">patient-hap1\n{sequence}\n", encoding="utf-8")
    metadata = {"schema_version": "1.0", "sequences": [{"record_schema_version": "1.0", "sequence_id": "patient-hap1",
        "sequence_role": "PATIENT_HAPLOTYPE", "display_label": "Patient sequence 1", "source_fasta_record_id": "SOURCE_HAP1",
        "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(), "sequence_length": 4, "original_orientation": "AS_STORED",
        "reference_alignment_strand": "NOT_ALIGNED", "display_orientation": "UNRESOLVED", "reverse_complement_required": None,
        "source_coordinates": None, "mapping_status": "NOT_ALIGNED", "warnings": []}]}
    (patient / "patient-sequences.metadata.json").write_text(json.dumps(metadata), encoding="utf-8")
    (patient / "assembly_record_mappings.json").write_text(json.dumps({"schema_version": "1.0", "coordinate_convention": "zero_based_half_open", "mappings": []}), encoding="utf-8")
    locus = tmp_path / "locus config.yaml"; locus.write_text('schema_version: "1.0"\n', encoding="utf-8")
    run_config = tmp_path / "run config.yaml"; run_config.write_text('schema_version: "1.0"\n', encoding="utf-8")
    identity = {"case_id": "C", "subject_id": "P", "sample_id": "S", "locus_id": "L"}
    run_normalization_stage({"run": identity, "locus_config": str(locus)}, run)
    package = tmp_path / "portable package"
    config = {"run": identity, "locus_config": str(locus), "_config_path": str(run_config),
              "case_package": {"package_root": str(package), "include_native_outputs": False,
                               "include_prepared_mini_bam": False, "include_alignment_artifacts": False,
                               "include_command_records": False}}
    before = {path.relative_to(run).as_posix(): (path.stat().st_size, sha256_file(path)) for path in run.rglob("*") if path.is_file()}
    build_case_package(config, run)
    after = {path.relative_to(run).as_posix(): (path.stat().st_size, sha256_file(path)) for path in run.rglob("*") if path.is_file()}
    return run, package, before, after


def test_builder_does_not_mutate_source_run(tmp_path):
    run, package, before, after = _built_package(tmp_path)
    assert before == after
    assert not (run / "00_manifest/pipeline-run-manifest.json").exists()
    assert validate_case_package(package)["valid"]


def test_moved_package_validates_after_original_run_is_deleted(tmp_path):
    run, package, _, _ = _built_package(tmp_path)
    moved = tmp_path / "different machine" / "case"
    moved.parent.mkdir(); shutil.move(package, moved); shutil.rmtree(run)
    report = validate_case_package(moved)
    assert report["valid"], report["issues"]


@pytest.mark.parametrize("relative", ["evidence/evidence-summary.json", "evidence/source-registry.json"])
def test_malformed_packaged_contract_is_detected(tmp_path, relative):
    _, package, _, _ = _built_package(tmp_path)
    document = json.loads((package / relative).read_text())
    document["unexpected_schema_field"] = True
    (package / relative).write_text(json.dumps(document), encoding="utf-8")
    report = validate_case_package(package)
    assert any(item["code"] == "SCHEMA_VALIDATION_ERROR" and item["relative_path"] == relative for item in report["issues"])


@pytest.mark.parametrize("kind,code", [("unknown", "UNKNOWN_MANIFEST_FILE_ID"), ("wrong-role", "WRONG_MANIFEST_FILE_ROLE")])
def test_manifest_reference_integrity(tmp_path, kind, code):
    _, package, _, _ = _built_package(tmp_path)
    path = package / "case-manifest.json"; manifest = json.loads(path.read_text())
    manifest["patient"]["fasta_file_id"] = ("missing-file" if kind == "unknown" else manifest["patient"]["metadata_file_id"])
    path.write_text(json.dumps(manifest), encoding="utf-8")
    assert any(item["code"] == code for item in validate_case_package(package)["issues"])
