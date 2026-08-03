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
        "VALIDATED_INPUTS": (root / "01_inputs" / "validated-inputs.json",),
        "PATIENT_SEQUENCE_METADATA": (root / "01_inputs" / "patient-sequences.metadata.json",),
        "PREPARED_MINI_BAM": (root / "02_prepared_bam" / "prepared.mini.bam",),
        "PREPARED_MINI_BAM_INDEX": (root / "02_prepared_bam" / "prepared.mini.bam.bai",),
        "BAM_VALIDATION_REPORT": (root / "02_prepared_bam" / "input_bam.validation.json",),
        "BAM_FLAGSTAT": (root / "02_prepared_bam" / "input_bam.flagstat.txt",),
        "BAM_IDXSTATS": (root / "02_prepared_bam" / "input_bam.idxstats.txt",),
        "ALIGNED_ASSEMBLY": (root / "03_assembly_alignment" / "assembly.aligned.sorted.bam",),
        "ALIGNED_ASSEMBLY_INDEX": (root / "03_assembly_alignment" / "assembly.aligned.sorted.bam.bai",),
        "ASSEMBLY_ALIGNMENT_SUMMARY": (root / "03_assembly_alignment" / "assembly_alignment.summary.json",),
        "ASSEMBLY_RECORD_MAPPINGS": (root / "03_assembly_alignment" / "assembly_record_mappings.json",),
        "PACKAGE_PATIENT_FASTA": (root / "03_assembly_alignment" / "patient-sequences.fasta",),
        "PACKAGE_PATIENT_METADATA": (root / "03_assembly_alignment" / "patient-sequences.metadata.json",),
        "VAMOS_READ_NATIVE_OUTPUTS": (root / "04_vamos_read" / "vamos-read.outputs.json",),
        "VAMOS_READ_RUN_METADATA": (root / "04_vamos_read" / "vamos-read.run.json",),
        "VAMOS_READ_NORMALIZED_EVIDENCE": (root / "04_vamos_read" / "vamos-read.normalized.json",),
        "VAMOS_CONTIG_NATIVE_OUTPUTS": (root / "05_vamos_contig" / "stage-outputs.json",),
        "VAMOS_CONTIG_RUN_METADATA": (root / "05_vamos_contig" / "stage-summary.json",),
        "VAMOS_CONTIG_NORMALIZED_EVIDENCE": (root / "05_vamos_contig" / "stage-normalized.json",),
        "STRAGLR_NATIVE_OUTPUTS": (root / "06_straglr" / "straglr.outputs.json",),
        "STRAGLR_RUN_METADATA": (root / "06_straglr" / "straglr.run.json",),
        "STRAGLR_NORMALIZED_EVIDENCE": (root / "06_straglr" / "straglr.normalized.json",),
        "TANDEM_GENOTYPES_ALIGNMENT_INPUTS": (root / "07_tandem_genotypes_preparation" / "alignment-inputs.json",),
        "TANDEM_GENOTYPES_ALIGNMENT_METADATA": (root / "07_tandem_genotypes_preparation" / "alignment-metadata.json",),
        "TANDEM_GENOTYPES_PREPARATION_SUMMARY": (root / "07_tandem_genotypes_preparation" / "preparation-summary.json",),
        "TANDEM_GENOTYPES_NATIVE_OUTPUTS": (root / "08_tandem_genotypes" / "stage-outputs.json",),
        "TANDEM_GENOTYPES_RUN_METADATA": (root / "08_tandem_genotypes" / "stage-summary.json",),
        "TANDEM_GENOTYPES_NORMALIZED_EVIDENCE": (root / "08_tandem_genotypes" / "stage-normalized.json",),
        "NATIVE_CALLER_OUTPUTS": (
            root / "04_vamos_read" / "vamos_read.vcf.gz",
            root / "05_vamos_contig" / "vamos_contig.vcf.gz",
            root / "06_straglr" / "straglr.tsv",
            root / "07_tandem_genotypes" / "tandem_genotypes.txt",
        ),
        "UNIFIED_NORMALIZED_EVIDENCE": (root / "09_normalized_evidence" / "normalized-evidence.json",),
        "UNIFIED_EVIDENCE_SUMMARY": (root / "09_normalized_evidence" / "evidence-summary.json",),
        "UNIFIED_SOURCE_REGISTRY": (root / "09_normalized_evidence" / "source-registry.json",),
        "UNIFIED_VALIDATION_REPORT": (root / "09_normalized_evidence" / "validation-report.json",),
        "NORMALIZED_EVIDENCE": (root / "09_normalized_evidence" / "normalized-evidence.json",),
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
                if tool_id in {"LASTDB", "LASTAL"}:
                    tool = detect_version(tool, pattern=r"(?i)(?:last(?:db|al)?[^0-9]*)?([0-9]{3,5})")
                elif tool_id == "TANDEM_GENOTYPES":
                    tool = detect_version(tool, pattern=r"(?i)(?:tandem-genotypes[^0-9]*)?v?([0-9]+(?:\.[0-9]+)+)")
                else:
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


def _native_tool(config: dict[str, Any], config_directory: Path, tool_id: str) -> Tool:
    executable, required = DEFAULT_TOOLS[tool_id]
    settings = config.get("tools", {}).get(TOOL_CONFIG_KEYS[tool_id], {})
    tool = resolve_tool(Tool(ToolId(tool_id), tool_id.replace("_", " ").title(),
        settings.get("executable", executable), settings.get("required", required)), config_directory)
    if not tool.resolved_executable:
        raise ConfigurationError(f"required tool {tool_id} is unavailable: {tool.configured_executable}")
    return detect_version(tool)


def _execute_implemented_stage(stage: StageDefinition, config: dict[str, Any], root: Path, config_directory: Path, overwrite: bool) -> None:
    if stage.stage_id == "00_validate_inputs":
        from .inputs import validate_inputs
        validate_inputs(config, root / "01_inputs")
    elif stage.stage_id == "01_prepare_bam":
        from .bam import prepare_bam
        prepare_bam(config["inputs"]["mini_bam"], config["inputs"]["mini_bam_index"], root / "02_prepared_bam",
            _native_tool(config, config_directory, "SAMTOOLS"), overwrite=overwrite)
    elif stage.stage_id == "02_align_assembly":
        from .alignment import align_assembly
        align_assembly(config, root / "03_assembly_alignment", _native_tool(config, config_directory, "MINIMAP2"),
            _native_tool(config, config_directory, "SAMTOOLS"), overwrite=overwrite)
    elif stage.stage_id in {"03_run_vamos_read", "04_run_vamos_contig"}:
        from .callers.vamos import run_vamos_stage
        run_vamos_stage(stage.stage_id, config, root, config_directory, overwrite=overwrite)
    elif stage.stage_id == "05_run_straglr":
        from .callers.straglr import run_straglr_stage
        run_straglr_stage(config, root, config_directory, overwrite=overwrite)
    elif stage.stage_id == "06_prepare_tandem_genotypes":
        from .last_alignment import prepare_tandem_genotypes_stage
        prepare_tandem_genotypes_stage(config, root, config_directory, overwrite=overwrite)
    elif stage.stage_id == "07_run_tandem_genotypes":
        from .callers.tandem_genotypes import run_tandem_genotypes_stage
        run_tandem_genotypes_stage(config, root, config_directory, overwrite=overwrite)
    elif stage.stage_id == "08_normalize_outputs":
        from .normalization import run_normalization_stage
        run_normalization_stage(config, root, overwrite=overwrite)


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
        if stage.stage_id == "05_run_straglr":
            # The catalog is an explicit locus resource rather than a filename-
            # inferred role. Its bytes therefore participate directly in resume.
            from .callers.straglr import _catalog
            catalog = _catalog(config)
            if catalog is not None:
                input_paths = (*input_paths, catalog)
        if stage.stage_id in {"06_prepare_tandem_genotypes", "07_run_tandem_genotypes"}:
            from .last_alignment import _resource
            repeat = _resource(config)
            if repeat is not None:
                input_paths = (*input_paths, repeat)
        if stage.stage_id == "07_run_tandem_genotypes":
            registry_path = root / "07_tandem_genotypes_preparation" / "alignment-inputs.json"
            if registry_path.is_file():
                registry = json.loads(registry_path.read_text(encoding="utf-8"))
                input_paths = (*input_paths, *(Path(item["path"]) for item in registry.get("alignments", [])))
        if stage.stage_id == "08_normalize_outputs":
            for registry_path in (root / "04_vamos_read" / "vamos-read.outputs.json", root / "05_vamos_contig" / "stage-outputs.json",
                                  root / "06_straglr" / "straglr.outputs.json", root / "08_tandem_genotypes" / "stage-outputs.json"):
                if registry_path.is_file():
                    registry = json.loads(registry_path.read_text(encoding="utf-8"))
                    input_paths = (*input_paths, *(Path(item["path"]) if Path(item["path"]).is_absolute() else (registry_path.parent / item["path"]).resolve() for item in registry.get("outputs", []) if item.get("path")))
        output_paths = tuple(path for role in stage.expected_output_roles for path in role_paths.get(role, ()))
        inputs = _identities(input_paths)
        prior_output_paths = tuple(Path(str(item["path"])) for item in (prior or {}).get("output_file_identities", []))
        outputs = _identities((*output_paths, *prior_output_paths))
        tools = _tool_identities(stage, config, config_directory, execution_mode)
        evaluated_reason = resume_eligibility(prior, digest, inputs, outputs, tools)
        if evaluated_reason is ResumeReason.RESUME_ALLOWED and any(not path.is_file() for path in input_paths):
            evaluated_reason = ResumeReason.INPUT_CHANGED

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
        implemented = stage.order <= 8
        status = StageStatus.DRY_RUN.value if dry_run else StageStatus.PLANNED.value
        warnings = (["Dry run: biological outputs and external commands were not created."] if dry_run and implemented else
            (["Caller-specific execution is deferred; this record describes scaffold planning only."] if not implemented else []))
        record = {
            "record_schema_version": "1.0",
            "stage_id": stage.stage_id,
            "status": status,
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
            "warnings": warnings,
            "failure": None,
            "resume_eligibility": {"eligible": False, "reason": reason.value},
        }
        atomic_write_json(record_path, record)
        if implemented and dry_run and stage.order in (3, 4, 5, 6, 7):
            if stage.order == 5:
                from .callers.straglr import run_straglr_stage
                run_straglr_stage(config, root, config_directory, overwrite=overwrite, dry_run=True)
            elif stage.order in (3, 4):
                from .callers.vamos import run_vamos_stage
                run_vamos_stage(stage.stage_id, config, root, config_directory, overwrite=overwrite, dry_run=True)
            elif stage.order == 6:
                from .last_alignment import prepare_tandem_genotypes_stage
                prepare_tandem_genotypes_stage(config, root, config_directory, overwrite=overwrite, dry_run=True)
            else:
                from .callers.tandem_genotypes import run_tandem_genotypes_stage
                run_tandem_genotypes_stage(config, root, config_directory, overwrite=overwrite, dry_run=True)
        if implemented and not dry_run:
            record["status"] = StageStatus.RUNNING.value
            record["completed_utc"] = None
            atomic_write_json(record_path, record)
            try:
                _execute_implemented_stage(stage, config, root, config_directory, overwrite)
                missing = [str(path) for path in output_paths if not path.is_file() or path.stat().st_size == 0]
                if missing:
                    raise RuntimeError(f"stage required outputs are missing or empty: {', '.join(missing)}")
                record["status"] = StageStatus.SUCCEEDED.value
                completed_output_paths = output_paths
                if stage.stage_id == "05_run_straglr":
                    registry = json.loads((root / "06_straglr" / "straglr.outputs.json").read_text(encoding="utf-8"))
                    completed_output_paths = (*output_paths, *(Path(item["path"]) for item in registry["outputs"]))
                elif stage.stage_id == "06_prepare_tandem_genotypes":
                    registry=json.loads((root/"07_tandem_genotypes_preparation"/"alignment-inputs.json").read_text())
                    metadata=json.loads((root/"07_tandem_genotypes_preparation"/"alignment-metadata.json").read_text())
                    completed_output_paths=(*output_paths,*(Path(item["path"]) for item in registry["alignments"]),*(Path(db["path"]) for item in metadata["records"] for db in item.get("database_file_ids",[])))
                elif stage.stage_id == "07_run_tandem_genotypes":
                    registry=json.loads((root/"08_tandem_genotypes"/"stage-outputs.json").read_text())
                    completed_output_paths=(*output_paths,*(Path(item["path"]) for item in registry["outputs"]))
                record["output_file_identities"] = _identities(completed_output_paths)
                execution_root = {1: root / "02_prepared_bam", 2: root / "03_assembly_alignment", 3: root / "04_vamos_read", 4: root / "05_vamos_contig", 5: root / "06_straglr", 6:root/"07_tandem_genotypes_preparation",7:root/"08_tandem_genotypes",8:root/"09_normalized_evidence"}.get(stage.order)
                record["command_record_paths"] = [str(path) for path in sorted(execution_root.glob("**/execution-records/*.json"))] if execution_root else []
                record["completed_utc"] = utc_now()
            except Exception as exc:
                record["status"] = StageStatus.FAILED.value
                record["completed_utc"] = utc_now()
                record["failure"] = {"error_type": type(exc).__name__, "message": str(exc)}
                atomic_write_json(record_path, record)
                raise
            atomic_write_json(record_path, record)
    return root
