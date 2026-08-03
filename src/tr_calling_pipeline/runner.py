"""Truthful scaffold planning with tool, option, and resume provenance."""

from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path
import sys
from typing import Any, Iterable

from .config import load_config
from .errors import ConfigurationError, StageResumeError
from .provenance import atomic_write_json, canonical_digest, file_identity, sha256_file, utc_now
from .stages import ResumeReason, StageDefinition, StageStatus, resume_eligibility, select_stages
from .tools import ExecutionMode, Tool, ToolId, detect_version, resolve_tool
from .version import pipeline_identity


DEFAULT_TOOLS = {
    "SAMTOOLS": ("samtools", True),
    "MINIMAP2": ("minimap2", True),
    "VAMOS": ("vamos", False),
    "STRAGLR": ("straglr.py", False),
    "LASTDB": ("lastdb", False),
    "LASTAL": ("lastal", False),
    "TANDEM_GENOTYPES": ("tandem-genotypes", False),
}
TOOL_CONFIG_KEYS = {tool_id: tool_id.lower() for tool_id in DEFAULT_TOOLS}


def _run_root(config: dict[str, Any]) -> Path:
    return Path(config["run"]["output_root"]) / f'{config["run"]["sample_id"]}_{config["run"]["locus_id"]}'


def _role_paths(config: dict[str, Any], root: Path) -> dict[str, tuple[Path, ...]]:
    inputs = config["inputs"]
    return {
        "ASSEMBLY_FASTA": (Path(inputs["assembly_fasta"]),),
        "MINI_BAM": (Path(inputs["mini_bam"]),),
        "MINI_BAM_INDEX": (Path(inputs["mini_bam_index"]),),
        "REFERENCE_FASTA": (Path(inputs["reference_fasta"]),),
        "REFERENCE_FASTA_INDEX": (Path(inputs["reference_fasta_index"]),),
        "LOCUS_CONFIG": (Path(config["locus_config"]),),
        "VALIDATED_INPUTS": (root / "00_manifest" / "validated-inputs.json",),
        "PREPARED_MINI_BAM": (root / "02_prepared_bam" / "original.mini.bam",),
        "PREPARED_MINI_BAM_INDEX": (root / "02_prepared_bam" / "original.mini.bam.bai",),
        "ALIGNED_ASSEMBLY": (root / "03_assembly_alignment" / "assembly.bam",),
        "ALIGNED_ASSEMBLY_INDEX": (root / "03_assembly_alignment" / "assembly.bam.bai",),
        "VAMOS_READ_NATIVE_OUTPUT": (root / "04_vamos_read" / "vamos_read.vcf.gz",),
        "VAMOS_CONTIG_NATIVE_OUTPUT": (root / "05_vamos_contig" / "vamos_contig.vcf.gz",),
        "STRAGLR_NATIVE_OUTPUT": (root / "06_straglr" / "straglr.tsv",),
        "TANDEM_GENOTYPES_ALIGNMENT_INPUT": (root / "07_tandem_genotypes" / "alignment.maf",),
        "TANDEM_GENOTYPES_NATIVE_OUTPUT": (root / "07_tandem_genotypes" / "tandem_genotypes.txt",),
        "NATIVE_CALLER_OUTPUTS": (
            root / "04_vamos_read" / "vamos_read.vcf.gz",
            root / "05_vamos_contig" / "vamos_contig.vcf.gz",
            root / "06_straglr" / "straglr.tsv",
            root / "07_tandem_genotypes" / "tandem_genotypes.txt",
        ),
        "NORMALIZED_EVIDENCE": (root / "08_normalized" / "caller_results.json",),
        "CASE_PACKAGE": (Path(config["case_package"]["package_root"]) / "case-manifest.json",),
        "VALIDATED_CASE_PACKAGE": (Path(config["case_package"]["package_root"]) / "case-manifest.json",),
    }


def _identities(paths: Iterable[Path]) -> list[dict[str, object]]:
    return [asdict(file_identity(path)) for path in paths if path.is_file()]


def _configured_mode(tool_config: dict[str, Any], override: str | None) -> ExecutionMode:
    return ExecutionMode(override or tool_config.get("execution_mode", ExecutionMode.NATIVE.value))


def _validate_container(config: dict[str, Any], config_directory: Path) -> tuple[str, str]:
    container = config.get("container") or {}
    apptainer = container.get("apptainer") or {}
    if not apptainer.get("executable") or not container.get("image"):
        raise ConfigurationError("APPTAINER execution requires container.apptainer.executable and container.image")
    image = Path(container["image"]).expanduser()
    image = image if image.is_absolute() else (config_directory / image).resolve()
    if not image.is_file():
        raise ConfigurationError(f"APPTAINER container image does not exist: {image}")
    return str(image), sha256_file(image)


def _tool_identities(
    stage: StageDefinition,
    config: dict[str, Any],
    config_directory: Path,
    execution_mode: str | None,
) -> list[dict[str, object]]:
    identities: list[dict[str, object]] = []
    pipeline = pipeline_identity()
    for tool_id in (*stage.required_tools, *stage.optional_tools):
        if tool_id == "PIPELINE":
            identities.append({"tool_id": "PIPELINE", "configured_executable": None, "resolved_executable": None, "detected_version": pipeline.pipeline_version, "execution_mode": ExecutionMode.NATIVE.value, "container_digest": None, "status": "AVAILABLE"})
            continue
        if tool_id == "PYTHON":
            identities.append({"tool_id": "PYTHON", "configured_executable": sys.executable, "resolved_executable": str(Path(sys.executable).resolve()), "detected_version": sys.version.split()[0], "execution_mode": ExecutionMode.NATIVE.value, "container_digest": None, "status": "AVAILABLE"})
            continue
        default_executable, default_required = DEFAULT_TOOLS[tool_id]
        tool_config = config.get("tools", {}).get(TOOL_CONFIG_KEYS[tool_id], {})
        mode = _configured_mode(tool_config, execution_mode)
        container_digest = None
        container_image = None
        if mode is ExecutionMode.APPTAINER:
            container_image, container_digest = _validate_container(config, config_directory)
        tool = Tool(
            tool_id=ToolId(tool_id),
            display_name=tool_id.replace("_", " ").title(),
            configured_executable=tool_config.get("executable", default_executable),
            required=tool_config.get("required", default_required),
            execution_mode=mode,
        )
        if mode is ExecutionMode.NATIVE:
            tool = resolve_tool(tool, config_directory)
            if tool.resolved_executable:
                tool = detect_version(tool)
            identity = tool.to_dict()
        else:
            identity = tool.to_dict()
            identity.update({
                "internal_executable": tool.configured_executable,
                "apptainer_executable": config["container"]["apptainer"]["executable"],
                "container_image": container_image,
            })
        identity["container_digest"] = container_digest
        identities.append(identity)
    return identities


def _stage_digest(stage: StageDefinition, config: dict[str, Any], mode: str | None) -> str:
    relevant_tools = {
        TOOL_CONFIG_KEYS[tool_id]: config.get("tools", {}).get(TOOL_CONFIG_KEYS[tool_id], {})
        for tool_id in (*stage.required_tools, *stage.optional_tools)
        if tool_id in TOOL_CONFIG_KEYS
    }
    return canonical_digest({
        "stage_id": stage.stage_id,
        "run": config["run"],
        "execution": config["execution"],
        "inputs": config["inputs"],
        "locus_config": config["locus_config"],
        "case_package": config["case_package"] if stage.order >= 9 else None,
        "tools": relevant_tools,
        "container": config.get("container", {}) if mode == ExecutionMode.APPTAINER.value else None,
        "execution_mode_override": mode,
    })


def _archive_invalidated(path: Path, prior: dict[str, Any], reason: ResumeReason) -> None:
    invalidated = dict(prior)
    invalidated["status"] = StageStatus.INVALIDATED.value
    invalidated["resume_eligibility"] = {"eligible": False, "reason": reason.value}
    invalidated["invalidated_utc"] = utc_now()
    archive = path.with_name(f"{path.stem}.invalidated.{invalidated['invalidated_utc'].replace(':', '')}.json")
    atomic_write_json(archive, invalidated)


def run(config_path, *, dry_run=False, start_stage=None, stop_stage=None, resume=True, overwrite=False, execution_mode=None):
    config = load_config(config_path)
    overwrite = bool(overwrite or config["run"].get("overwrite", False))
    config_directory = Path(config["_config_path"]).parent
    stages = select_stages(start_stage, stop_stage)
    root = _run_root(config)
    records = root / "00_manifest" / "stages"
    records.mkdir(parents=True, exist_ok=True)
    role_paths = _role_paths(config, root)

    # Validate a global override even when the selected stage has no external tool.
    if execution_mode == ExecutionMode.APPTAINER.value:
        _validate_container(config, config_directory)

    for stage in stages:
        record_path = records / f"{stage.stage_id}.json"
        prior = json.loads(record_path.read_text(encoding="utf-8")) if record_path.is_file() else None
        digest = _stage_digest(stage, config, execution_mode)
        input_paths = tuple(path for role in stage.required_input_roles for path in role_paths.get(role, ()))
        output_paths = tuple(path for role in stage.expected_output_roles for path in role_paths.get(role, ()))
        inputs = _identities(input_paths)
        outputs = _identities(output_paths)
        tools = _tool_identities(stage, config, config_directory, execution_mode)
        evaluated_reason = resume_eligibility(prior, digest, inputs, outputs, tools)
        if evaluated_reason is ResumeReason.RESUME_ALLOWED and any(not path.is_file() for path in input_paths):
            evaluated_reason = ResumeReason.INPUT_CHANGED
        if evaluated_reason is ResumeReason.RESUME_ALLOWED and any(not path.is_file() for path in output_paths):
            evaluated_reason = ResumeReason.OUTPUT_MISSING

        if resume and evaluated_reason is ResumeReason.RESUME_ALLOWED:
            now = utc_now()
            skip_record = {
                **prior,
                "status": StageStatus.SKIPPED.value,
                "started_utc": now,
                "completed_utc": now,
                "duration_seconds": 0.0,
                "resume_eligibility": {"eligible": True, "reason": evaluated_reason.value},
                "overwrite": overwrite,
            }
            # Preserve the canonical successful record so a later run can resume again.
            atomic_write_json(records / f"{stage.stage_id}.skip.json", skip_record)
            continue

        reason = evaluated_reason if resume else ResumeReason.RESUME_DISABLED
        if prior and not resume and not overwrite:
            raise StageResumeError(
                f"stage record already exists for {stage.stage_id}; use --resume or --overwrite",
                stage_id=stage.stage_id,
            )
        if prior and resume:
            _archive_invalidated(record_path, prior, reason)

        now = utc_now()
        effective_modes = sorted({str(tool["execution_mode"]) for tool in tools})
        record = {
            "record_schema_version": "1.0",
            "stage_id": stage.stage_id,
            "status": StageStatus.DRY_RUN.value if dry_run else StageStatus.PLANNED.value,
            "started_utc": now,
            "completed_utc": now,
            "duration_seconds": 0.0,
            "configuration_digest": digest,
            "input_file_identities": inputs,
            "output_file_identities": outputs,
            "tool_identities": tools,
            "execution_mode": execution_mode or (effective_modes[0] if len(effective_modes) == 1 else "NATIVE"),
            "overwrite": overwrite,
            "command_record_paths": [],
            "warnings": ["Caller-specific execution is deferred; this record describes scaffold planning only."],
            "failure": None,
            "resume_eligibility": {"eligible": False, "reason": reason.value},
        }
        atomic_write_json(record_path, record)
    return root
