import json
from pathlib import Path

import pytest

from tr_calling_pipeline.case_package import _safe_relative, reject_secret_bearing_config
from tr_calling_pipeline.case_package_validation import validate_case_package
from tr_calling_pipeline.config import ConfigurationError


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
    assert any(issue["code"] == "INVALID_MANIFEST" for issue in report["issues"])


def test_unexpected_file_is_an_error(tmp_path):
    # A deliberately incomplete package still demonstrates the strict extra-file policy.
    (tmp_path / "extra.txt").write_text("unexpected", encoding="utf-8")
    report = validate_case_package(tmp_path)
    assert any(issue["code"] == "UNEXPECTED_UNREGISTERED_FILE" for issue in report["issues"])
