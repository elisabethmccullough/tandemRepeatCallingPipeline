from copy import deepcopy
from dataclasses import asdict
import json
from pathlib import Path
import stat
import sys

import pytest
import yaml

from tr_calling_pipeline.cli import main
from tr_calling_pipeline.errors import ConfigurationError, StageResumeError
from tr_calling_pipeline.provenance import file_identity
from tr_calling_pipeline.runner import run
from tr_calling_pipeline.stages import STAGES
BASE = yaml.safe_load((Path(__file__).parents[1] / "config/example.yaml").read_text())


def make_run(tmp_path: Path) -> tuple[Path, Path, Path, Path]:
    document = deepcopy(BASE)
    document["run"]["output_root"] = "results"
    document["case_package"]["package_root"] = "package"
    for key in ("assembly_fasta", "mini_bam", "mini_bam_index", "reference_fasta", "reference_fasta_index"):
        path = tmp_path / "inputs" / key
        path.parent.mkdir(exist_ok=True)
        path.write_text(f"original {key}\n")
        document["inputs"][key] = str(path.relative_to(tmp_path))
    locus = tmp_path / "locus.yaml"
    locus.write_text("synthetic locus identity\n")
    document["locus_config"] = locus.name
    tool = tmp_path / "fake samtools.py"
    tool.write_text(f"#!{sys.executable}\nprint('samtools 1.0')\n")
    tool.chmod(tool.stat().st_mode | stat.S_IXUSR)
    document["tools"]["samtools"]["executable"] = f"./{tool.name}"
    config = tmp_path / "run config.yaml"
    config.write_text(yaml.safe_dump(document))
    root = tmp_path / "results" / "HG00438_HTT"
    output = root / "02_prepared_bam" / "original.mini.bam"
    return config, root, output, tool


def successful_prior(config: Path, root: Path, output: Path) -> Path:
    run(config, dry_run=True, start_stage="01_prepare_bam", stop_stage="01_prepare_bam")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("prepared bam\n")
    output_index = output.with_suffix(".bam.bai")
    output_index.write_text("prepared bam index\n")
    record_path = root / "00_manifest" / "stages" / "01_prepare_bam.json"
    record = json.loads(record_path.read_text())
    record["status"] = "SUCCEEDED"
    record["output_file_identities"] = [asdict(file_identity(output)), asdict(file_identity(output_index))]
    record_path.write_text(json.dumps(record, indent=2) + "\n")
    return record_path


def test_resume_skips_and_no_resume_requires_overwrite(tmp_path):
    config, root, output, _ = make_run(tmp_path)
    record_path = successful_prior(config, root, output)
    run(config, dry_run=True, start_stage="01_prepare_bam", stop_stage="01_prepare_bam", resume=True)
    assert json.loads(record_path.read_text())["status"] == "SUCCEEDED"
    skip = json.loads(record_path.with_name("01_prepare_bam.skip.json").read_text())
    assert skip["status"] == "SKIPPED"
    assert skip["resume_eligibility"] == {"eligible": True, "reason": "RESUME_ALLOWED"}
    with pytest.raises(StageResumeError, match="--overwrite"):
        run(config, dry_run=True, start_stage="01_prepare_bam", stop_stage="01_prepare_bam", resume=False)
    run(config, dry_run=True, start_stage="01_prepare_bam", stop_stage="01_prepare_bam", resume=False, overwrite=True)
    replaced = json.loads(record_path.read_text())
    assert replaced["status"] == "DRY_RUN"
    assert replaced["overwrite"] is True
    assert replaced["resume_eligibility"]["reason"] == "RESUME_DISABLED"


@pytest.mark.parametrize("changed,reason", [("input", "INPUT_CHANGED"), ("output", "OUTPUT_CHANGED"), ("configuration", "CONFIGURATION_CHANGED"), ("version", "TOOL_VERSION_CHANGED")])
def test_resume_rejection_reasons(tmp_path, changed, reason):
    config, root, output, tool = make_run(tmp_path)
    record_path = successful_prior(config, root, output)
    if changed == "input":
        (tmp_path / "inputs" / "mini_bam").write_text("corrupt input\n")
    elif changed == "output":
        output.write_text("corrupt output\n")
    elif changed == "configuration":
        document = yaml.safe_load(config.read_text())
        document["execution"]["threads"] += 1
        config.write_text(yaml.safe_dump(document))
    else:
        tool.write_text(f"#!{sys.executable}\nprint('samtools 2.0')\n")
        tool.chmod(tool.stat().st_mode | stat.S_IXUSR)
    run(config, dry_run=True, start_stage="01_prepare_bam", stop_stage="01_prepare_bam", resume=True)
    record = json.loads(record_path.read_text())
    assert record["status"] == "DRY_RUN"
    assert record["resume_eligibility"]["reason"] == reason
    assert list(record_path.parent.glob("01_prepare_bam.invalidated.*.json"))


def test_execution_mode_and_cli_options_are_observable(tmp_path):
    config, root, _, _ = make_run(tmp_path)
    assert main(["run", "--config", str(config), "--dry-run", "--no-resume", "--overwrite", "--execution-mode", "NATIVE", "--start-stage", "01_prepare_bam", "--stop-stage", "01_prepare_bam"]) == 0
    record = json.loads((root / "00_manifest/stages/01_prepare_bam.json").read_text())
    assert record["execution_mode"] == "NATIVE"
    assert record["overwrite"] is True
    assert record["resume_eligibility"]["reason"] == "RESUME_DISABLED"
    assert all(tool["execution_mode"] == "NATIVE" for tool in record["tool_identities"])


def test_apptainer_override_requires_complete_configuration(tmp_path):
    config, _, _, _ = make_run(tmp_path)
    with pytest.raises(ConfigurationError, match="container.image"):
        run(config, dry_run=True, execution_mode="APPTAINER", start_stage="01_prepare_bam", stop_stage="01_prepare_bam")


def test_stage_registry_has_workflow_roles_and_tools():
    stages = {stage.stage_id: stage for stage in STAGES}
    assert stages["01_prepare_bam"].required_input_roles == ("MINI_BAM", "MINI_BAM_INDEX")
    assert stages["01_prepare_bam"].required_tools == ("SAMTOOLS",)
    assert stages["02_align_assembly"].required_tools == ("MINIMAP2", "SAMTOOLS")
    assert stages["05_run_straglr"].optional_tools == ("STRAGLR",)
    assert stages["08_normalize_outputs"].expected_output_roles == (
        "UNIFIED_NORMALIZED_EVIDENCE", "UNIFIED_EVIDENCE_SUMMARY",
        "UNIFIED_SOURCE_REGISTRY", "UNIFIED_VALIDATION_REPORT")
    assert stages["09_build_case_package"].expected_output_roles == ("CASE_PACKAGE",)
    assert stages["10_validate_case_package"].required_input_roles == ("CASE_PACKAGE",)
