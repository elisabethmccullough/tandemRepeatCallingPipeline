import json
from pathlib import Path

from tr_calling_pipeline.readiness import check_release_readiness
from tr_calling_pipeline.schema_validation import validate
from tr_calling_pipeline.verification import VERIFICATION_LEVELS, schema_directory, validate_fixtures, validate_schemas
from tr_calling_pipeline.synthetic_demo import _prepare_workspace, run_demo, run_package_portability_demo
from tr_calling_pipeline.provenance import sha256_file
from tr_calling_pipeline.runner import run


def test_verification_vocabulary_and_schema():
    assert "LABORATORY_VERIFIED" in VERIFICATION_LEVELS
    schema = json.loads((schema_directory() / "component-verification-status.schema.json").read_text())
    validate({"record_schema_version": "1.0", "component_id": "vamos-adapter", "component_type": "TOOL_ADAPTER",
              "implementation_status": "DEVELOPMENT_GATED", "verification_level": "SYNTHETIC_INTEGRATION_TESTED",
              "verified_versions": [], "test_evidence": ["tests/test_vamos.py"],
              "known_limitations": ["fake tool only"], "required_follow_up": ["real-tool smoke test"],
              "updated_utc": "2026-08-03T00:00:00Z"}, schema)


def test_all_packaged_schemas_parse_and_have_unique_ids():
    report = validate_schemas()
    assert report["schema_count"] >= 30
    assert report["valid"], report["errors"]
    assert validate_fixtures()["valid"]


def test_readiness_is_truthful(monkeypatch):
    status, limitations = check_release_readiness(Path(__file__).parents[1])
    assert status == "NOT_READY"
    assert any("hosted CI" in item for item in limitations)
    monkeypatch.setenv("TR_PIPELINE_FULL_DEMO_VERIFIED", "1")
    monkeypatch.setenv("TR_PIPELINE_WHEEL_VERIFIED", "1")
    monkeypatch.setenv("TR_PIPELINE_HOSTED_CI_VERIFIED", "1")
    status, limitations = check_release_readiness(Path(__file__).parents[1])
    assert status == "READY_WITH_LIMITATIONS"
    assert any("laboratory" in item for item in limitations)


def test_no_adapter_claims_laboratory_verification():
    roots = [Path(__file__).parents[1] / "docs", Path(__file__).parents[1] / "tests" / "fixtures"]
    claims = []
    for root in roots:
        for path in root.rglob("*.json"):
            if '"verification_level": "LABORATORY_VERIFIED"' in path.read_text(encoding="utf-8", errors="ignore"):
                claims.append(path)
    assert not claims


def test_package_portability_demo_moves_package_and_removes_source(tmp_path):
    report = run_package_portability_demo(tmp_path / "package path with spaces")
    assert report["valid_before_move"] and report["valid_after_move"] and report["valid_after_source_removal"]
    assert report["scope"] == "PACKAGE_PORTABILITY_ONLY"
    assert not (tmp_path / "package path with spaces" / "temporary package-only run").exists()


def test_full_synthetic_runner_demo(tmp_path):
    report = run_demo(tmp_path / "full runner path with spaces")
    assert report["stage_count"] == report["resumed_stage_count"] == 11
    assert all(count > 0 for count in report["caller_records_by_source"].values())
    assert report["valid_before_move"] and report["valid_after_move"] and report["valid_after_source_removal"]
    assert report["real_tools_executed"] is False
    package = Path(report["package"])
    evidence = json.loads((package / "evidence/normalized-evidence.json").read_text())
    read_records = [item for item in evidence["records"] if item["analysis_source"] == "RAW_READS"]
    assembly_records = [item for item in evidence["records"] if item["analysis_source"] == "ASSEMBLED_CONTIG"]
    assert all(item["assignment_state"] == "UNASSIGNED" and item["associated_sequence_id"] is None for item in read_records)
    assert {item["associated_sequence_id"] for item in assembly_records} == {"synthetic-sequence-1", "synthetic-sequence-2"}
    assert not any("consensus" in key.lower() or "classification" in key.lower() for item in evidence["records"] for key in item)
    tools = json.loads((package / "provenance/tools.json").read_text())["tools"]
    external = [item for item in tools if item["tool_id"] not in {"PIPELINE", "PYTHON"}]
    assert external and all(item["executable_kind"] == "FAKE_TEST_TOOL" for item in external)
    assert all(item["verification_level"] == "SYNTHETIC_INTEGRATION_TESTED" for item in external)
    command_records = [json.loads(path.read_text()) for path in (package / "provenance/commands").glob("*.json")]
    assert command_records and all(item["status"] == "SUCCEEDED" for item in command_records)
    assert any("fake-" in Path(item["argv"][0]).name for item in command_records)
    sources = json.loads((package / "evidence/source-registry.json").read_text())["sources"]
    manifest = json.loads((package / "case-manifest.json").read_text())
    package_ids = {item["source_file_id"]: item["file_id"] for item in manifest["native_sources"]}
    for source in sources:
        native = package / next(item["relative_path"] for item in manifest["files"] if item["file_id"] == package_ids[source["file_id"]])
        assert sha256_file(native) == source["sha256"]


def test_full_runner_upstream_change_invalidates_downstream(tmp_path):
    output = tmp_path / "invalidation path with spaces"; output.mkdir()
    config_path, package, _ = _prepare_workspace(output)
    root = Path(run(config_path))
    bam = config_path.parent / "fixtures/mini.bam"
    bam.write_text(bam.read_text() + "CHANGED\n")
    run(config_path, overwrite=True)
    archives = {path.name for path in (root / "00_manifest/stages").glob("*.invalidated.*.json")}
    assert any(name.startswith("01_prepare_bam") for name in archives)
    assert any(name.startswith("08_normalize_outputs") for name in archives)
    assert any(name.startswith("09_build_case_package") for name in archives)
    assert package.is_dir()
