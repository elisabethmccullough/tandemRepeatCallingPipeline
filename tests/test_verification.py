import json
from pathlib import Path

from tr_calling_pipeline.readiness import check_release_readiness
from tr_calling_pipeline.schema_validation import validate
from tr_calling_pipeline.verification import VERIFICATION_LEVELS, schema_directory, validate_schemas
from tr_calling_pipeline.synthetic_demo import run_demo


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


def test_readiness_is_truthful():
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


def test_synthetic_demo_moves_package_and_removes_source(tmp_path):
    report = run_demo(tmp_path / "path with spaces")
    assert report["valid_before_move"] and report["valid_after_move"] and report["valid_after_source_removal"]
    assert report["external_tools_executed"] is False
    assert not (tmp_path / "path with spaces" / "temporary synthetic run").exists()
