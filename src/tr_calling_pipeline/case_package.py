"""Construction of the portable, non-interpretive case-package contract.

The builder deliberately copies evidence byte-for-byte.  It inventories and
validates evidence; it never reconciles callers or derives a biological call.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import mimetypes
import os
from pathlib import Path, PurePosixPath
import shutil
import tempfile
from typing import Any, Iterable, Mapping
import uuid

import yaml

from .config import ConfigurationError
from .provenance import atomic_write_json, canonical_digest, sha256_file, utc_now
from .version import pipeline_identity


CONTRACT_VERSION = "1.0"
CONTROL_FILES = {"case-manifest.json", "checksums/sha256sums.json", "package-validation.json"}
SECRET_KEYS = {"token", "password", "secret", "api_key", "apikey", "credential", "credentials"}


@dataclass(frozen=True)
class CasePackageFile:
    file_id: str
    role: str
    relative_path: str
    media_type: str
    sha256: str
    size_bytes: int
    required: bool
    source_stage: str | None = None
    source_file_id: str | None = None
    source_path: str | None = None
    caller: str | None = None
    analysis_source: str | None = None
    associated_sequence_id: str | None = None


@dataclass(frozen=True)
class CasePackageValidationIssue:
    issue_id: str
    severity: str
    code: str
    message: str
    relative_path: str | None = None
    file_id: str | None = None
    record_id: str | None = None
    caller: str | None = None
    associated_sequence_id: str | None = None


@dataclass(frozen=True)
class CasePackageValidationReport:
    record_schema_version: str
    package_id: str | None
    valid: bool
    validated_utc: str
    error_count: int
    warning_count: int
    issues: tuple[CasePackageValidationIssue, ...]
    validated_files: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        value = asdict(self)
        value["issues"] = [asdict(issue) for issue in self.issues]
        value["validated_files"] = list(self.validated_files)
        return value


def _json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ConfigurationError(f"Cannot read JSON {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ConfigurationError(f"JSON root must be an object: {path}")
    return value


def _safe_relative(value: str) -> PurePosixPath:
    path = PurePosixPath(value)
    if not value or path.is_absolute() or ".." in path.parts or "\\" in value:
        raise ConfigurationError(f"unsafe package-relative path: {value}")
    return path


def _assert_regular(path: Path) -> None:
    if path.is_symlink() or not path.is_file():
        raise ConfigurationError(f"package sources must be regular files, not links or special files: {path}")


def _secret_paths(value: Any, prefix: str = "") -> list[str]:
    found: list[str] = []
    if isinstance(value, Mapping):
        for key, child in value.items():
            name = str(key).lower().replace("-", "_")
            location = f"{prefix}.{key}" if prefix else str(key)
            if name in SECRET_KEYS or any(part in SECRET_KEYS for part in name.split("_")):
                found.append(location)
            found.extend(_secret_paths(child, location))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            found.extend(_secret_paths(child, f"{prefix}[{index}]"))
    return found


def reject_secret_bearing_config(path: str | Path) -> None:
    source = Path(path)
    try:
        value = yaml.safe_load(source.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise ConfigurationError(f"cannot inspect configuration {source}: {exc}") from exc
    unsafe = sorted(set(_secret_paths(value)))
    if unsafe:
        raise ConfigurationError("configuration contains secret-like keys and will not be packaged: " + ", ".join(unsafe))


def _fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    current: str | None = None
    for number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            current = line[1:].split()[0]
            if not current or current in records:
                raise ConfigurationError(f"duplicate or empty FASTA record at line {number}")
            records[current] = []
        elif current is None:
            raise ConfigurationError(f"sequence before FASTA header at line {number}")
        else:
            records[current].append(line.upper())
    return {name: "".join(parts) for name, parts in records.items()}


def verify_patient_data(fasta_path: Path, metadata_path: Path) -> tuple[str, ...]:
    records = _fasta(fasta_path)
    metadata = _json(metadata_path)
    sequences = metadata.get("sequences")
    if not isinstance(sequences, list):
        raise ConfigurationError("patient metadata sequences must be an array")
    ids = [item.get("sequence_id") for item in sequences]
    if len(ids) != len(set(ids)):
        raise ConfigurationError("duplicate patient sequence ID")
    if set(ids) != set(records):
        raise ConfigurationError(f"patient FASTA/metadata sequence IDs differ: FASTA={sorted(records)}, metadata={sorted(ids)}")
    for item in sequences:
        sequence = records[item["sequence_id"]]
        digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()
        if item.get("sequence_sha256") != digest:
            raise ConfigurationError(f"biological sequence checksum mismatch: {item['sequence_id']}")
        if not item.get("source_fasta_record_id"):
            raise ConfigurationError(f"source FASTA record ID missing: {item['sequence_id']}")
    return tuple(ids)


def _file_id(role: str, relative_path: str, source_file_id: str | None = None) -> str:
    seed = source_file_id or f"{role}:{relative_path}"
    return "file-" + hashlib.sha256(seed.encode()).hexdigest()[:20]


def _native_destination(source: Mapping[str, Any]) -> str:
    caller = str(source.get("caller", "unknown")).upper()
    analysis = source.get("analysis_source")
    stem = {("VAMOS", "RAW_READS"): "vamos-read", ("VAMOS", "ASSEMBLED_CONTIG"): "vamos-contig",
            ("STRAGLR", "RAW_READS"): "straglr", ("TANDEM_GENOTYPES", "ASSEMBLED_CONTIG"): "tandem-genotypes"}.get((caller, analysis))
    if stem is None:
        raise ConfigurationError(f"invalid caller/source combination: {caller}/{analysis}")
    sequence = source.get("associated_sequence_id")
    directory = PurePosixPath("native", stem)
    if analysis == "ASSEMBLED_CONTIG":
        if not sequence:
            raise ConfigurationError(f"assembled-contig source lacks sequence association: {source.get('file_id')}")
        directory /= str(sequence)
    original = Path(str(source["path"])).name
    prefix = hashlib.sha256(str(source["file_id"]).encode()).hexdigest()[:10]
    return (directory / f"{prefix}-{original}").as_posix()


def _source_path(value: str, registry_path: Path, run_root: Path) -> Path:
    candidate = Path(value)
    if candidate.is_absolute():
        return candidate
    options = ((registry_path.parent / candidate), (run_root / candidate))
    return next((path.resolve() for path in options if path.is_file()), options[0].resolve())


def _copy(work: Path, source: Path, relative: str, role: str, *, required: bool = True,
          source_stage: str | None = None, source_file_id: str | None = None,
          caller: str | None = None, analysis_source: str | None = None,
          associated_sequence_id: str | None = None, media_type: str | None = None) -> CasePackageFile:
    _assert_regular(source)
    relative = _safe_relative(relative).as_posix()
    target = work.joinpath(*PurePosixPath(relative).parts)
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, target)
    if sha256_file(source) != sha256_file(target):
        raise ConfigurationError(f"copy verification failed: {source}")
    return CasePackageFile(_file_id(role, relative, source_file_id), role, relative,
        media_type or mimetypes.guess_type(source.name)[0] or "application/octet-stream",
        sha256_file(target), target.stat().st_size, required, source_stage, source_file_id,
        str(source), caller, analysis_source, associated_sequence_id)


def _write_summary(work: Path, relative: str, value: Any, role: str) -> CasePackageFile:
    target = work / relative
    atomic_write_json(target, value)
    return CasePackageFile(_file_id(role, relative), role, relative, "application/json", sha256_file(target),
                           target.stat().st_size, True, "09_build_case_package")


def plan_case_package(config: Mapping[str, Any], run_root: str | Path) -> dict[str, Any]:
    root = Path(run_root)
    options = {"include_prepared_mini_bam": True, "include_alignment_artifacts": True,
               "include_native_outputs": True, "include_command_records": True, **config.get("case_package", {})}
    required = {
        "PACKAGE_PATIENT_FASTA": root / "03_assembly_alignment/patient-sequences.fasta",
        "PACKAGE_PATIENT_METADATA": root / "03_assembly_alignment/patient-sequences.metadata.json",
        "ASSEMBLY_RECORD_MAPPINGS": root / "03_assembly_alignment/assembly_record_mappings.json",
        "UNIFIED_NORMALIZED_EVIDENCE": root / "09_normalized_evidence/normalized-evidence.json",
        "UNIFIED_EVIDENCE_SUMMARY": root / "09_normalized_evidence/evidence-summary.json",
        "UNIFIED_SOURCE_REGISTRY": root / "09_normalized_evidence/source-registry.json",
        "UNIFIED_VALIDATION_REPORT": root / "09_normalized_evidence/validation-report.json",
        "RUN_CONFIG": Path(str(config["_config_path"])), "LOCUS_CONFIG": Path(str(config["locus_config"])),
    }
    return {"status": "PLANNED", "package_root": str(config["case_package"]["package_root"]),
            "required": {key: str(path) for key, path in required.items()},
            "missing_required": sorted(key for key, path in required.items() if not path.is_file()), "options": options}


def build_case_package(config: Mapping[str, Any], run_root: str | Path, *, overwrite: bool = False,
                       dry_run: bool = False) -> dict[str, Any]:
    """Build, independently validate, and atomically publish one package."""
    plan = plan_case_package(config, run_root)
    if dry_run:
        return {**plan, "status": "DRY_RUN"}
    if plan["missing_required"]:
        raise ConfigurationError("missing required package inputs: " + ", ".join(plan["missing_required"]))
    required = {key: Path(value) for key, value in plan["required"].items()}
    reject_secret_bearing_config(required["RUN_CONFIG"])
    reject_secret_bearing_config(required["LOCUS_CONFIG"])
    sequence_ids = verify_patient_data(required["PACKAGE_PATIENT_FASTA"], required["PACKAGE_PATIENT_METADATA"])
    evidence = _json(required["UNIFIED_NORMALIZED_EVIDENCE"])
    report = _json(required["UNIFIED_VALIDATION_REPORT"])
    identity = config["run"]
    for field in ("case_id", "subject_id", "sample_id", "locus_id"):
        if evidence.get(field) != identity[field]:
            raise ConfigurationError(f"unified evidence {field} does not match run configuration")
    if report.get("valid") is not True:
        raise ConfigurationError("Stage 08 validation report is not valid")

    destination = Path(str(config["case_package"]["package_root"]))
    if destination.exists() and not overwrite:
        raise ConfigurationError(f"case package already exists: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix=f".{destination.name}.work-", dir=destination.parent))
    backup: Path | None = None
    try:
        files: list[CasePackageFile] = []
        fixed = (("PACKAGE_PATIENT_FASTA", "patient/patient-sequences.fasta", "PATIENT_FASTA"),
                 ("PACKAGE_PATIENT_METADATA", "patient/patient-sequences.metadata.json", "PATIENT_METADATA"),
                 ("ASSEMBLY_RECORD_MAPPINGS", "patient/assembly-record-mappings.json", "ASSEMBLY_RECORD_MAPPINGS"),
                 ("UNIFIED_NORMALIZED_EVIDENCE", "evidence/normalized-evidence.json", "NORMALIZED_EVIDENCE"),
                 ("UNIFIED_EVIDENCE_SUMMARY", "evidence/evidence-summary.json", "EVIDENCE_SUMMARY"),
                 ("UNIFIED_SOURCE_REGISTRY", "evidence/source-registry.json", "SOURCE_REGISTRY"),
                 ("UNIFIED_VALIDATION_REPORT", "evidence/normalization-validation-report.json", "NORMALIZATION_VALIDATION_REPORT"),
                 ("RUN_CONFIG", "config/run-config.yaml", "RUN_CONFIG"), ("LOCUS_CONFIG", "config/locus-config.yaml", "LOCUS_CONFIG"))
        for key, relative, role in fixed:
            files.append(_copy(work, required[key], relative, role, source_stage="08_normalize_outputs" if key.startswith("UNIFIED") else None))

        run_root_path = Path(run_root)
        registry = _json(required["UNIFIED_SOURCE_REGISTRY"])
        native_sources: list[dict[str, Any]] = []
        if plan["options"]["include_native_outputs"]:
            for source in sorted(registry.get("sources", []), key=lambda item: str(item.get("file_id"))):
                if source.get("associated_sequence_id") not in (None, *sequence_ids):
                    raise ConfigurationError(f"unknown native source sequence: {source.get('associated_sequence_id')}")
                origin = _source_path(str(source["path"]), required["UNIFIED_SOURCE_REGISTRY"], run_root_path)
                _assert_regular(origin)
                if sha256_file(origin) != source["sha256"] or origin.stat().st_size != source["size_bytes"]:
                    raise ConfigurationError(f"native source identity mismatch: {source['file_id']}")
                relative = _native_destination(source)
                entry = _copy(work, origin, relative, "CALLER_NATIVE_OUTPUT", required=False,
                    source_stage="08_normalize_outputs", source_file_id=source["file_id"], caller=source["caller"],
                    analysis_source=source["analysis_source"], associated_sequence_id=source.get("associated_sequence_id"),
                    media_type=source.get("media_type"))
                files.append(entry); native_sources.append({"file_id": entry.file_id, "source_file_id": source["file_id"]})

        options = plan["options"]
        if options["include_prepared_mini_bam"]:
            for origin, relative, role in ((run_root_path / "02_prepared_bam/prepared.mini.bam", "reads/prepared-mini.bam", "PREPARED_MINI_BAM"),
                                           (run_root_path / "02_prepared_bam/prepared.mini.bam.bai", "reads/prepared-mini.bam.bai", "PREPARED_MINI_BAM_INDEX")):
                if origin.is_file(): files.append(_copy(work, origin, relative, role, required=False, source_stage="01_prepare_bam"))
        if options["include_alignment_artifacts"]:
            for base, target in ((run_root_path / "03_assembly_alignment", "alignments/assembly"),
                                 (run_root_path / "07_tandem_genotypes_preparation", "alignments/last")):
                if base.is_dir():
                    for origin in sorted(path for path in base.rglob("*") if path.is_file() and not path.is_symlink()):
                        if origin in required.values(): continue
                        relative = f"{target}/{origin.relative_to(base).as_posix()}"
                        files.append(_copy(work, origin, relative, "ALIGNMENT_ARTIFACT", required=False))

        manifest_root = run_root_path / "00_manifest"
        pipeline_manifest = manifest_root / "pipeline-run-manifest.json"
        if pipeline_manifest.is_file():
            files.append(_copy(work, pipeline_manifest, "provenance/pipeline-run-manifest.json", "PIPELINE_MANIFEST"))
        else:
            # Focused construction must not repair or otherwise mutate its source
            # run.  Synthesize the missing provenance directly in the private
            # package work tree instead.
            files.append(_write_summary(work, "provenance/pipeline-run-manifest.json",
                {"record_schema_version": "1.0", "pipeline": asdict(pipeline_identity()), "run": identity,
                 "warning": "Package-local summary: the source run had no pipeline manifest."}, "PIPELINE_MANIFEST"))
        stages = [_json(path) for path in sorted((manifest_root / "stages").glob("*.json")) if ".invalidated." not in path.name]
        tools = sorted({json.dumps(tool, sort_keys=True): tool for stage in stages for tool in stage.get("tool_identities", [])}.values(), key=lambda x: str(x.get("tool_id")))
        files.append(_write_summary(work, "provenance/stages.json", {"record_schema_version": "1.0", "stages": stages}, "STAGE_PROVENANCE"))
        files.append(_write_summary(work, "provenance/tools.json", {"record_schema_version": "1.0", "tools": tools}, "TOOL_PROVENANCE"))
        files.append(_write_summary(work, "provenance/source-files.json", registry, "SOURCE_FILE_PROVENANCE"))
        if options["include_command_records"]:
            for origin in sorted(run_root_path.glob("**/execution-records/*.json")):
                relative = f"provenance/commands/{hashlib.sha256(str(origin).encode()).hexdigest()[:10]}-{origin.name}"
                files.append(_copy(work, origin, relative, "COMMAND_RECORD", required=False))

        files.sort(key=lambda item: item.relative_path)
        package_id = "case-" + canonical_digest({"identity": identity, "files": [(x.role, x.sha256, x.relative_path) for x in files],
                                                   "options": {k: options[k] for k in options if k.startswith("include_")}})[:24]
        by_role = {item.role: item.file_id for item in files}
        manifest = {"record_schema_version": "1.0", "package_contract_version": CONTRACT_VERSION,
            "package_id": package_id, **{key: identity[key] for key in ("case_id", "subject_id", "sample_id", "locus_id")},
            "created_utc": utc_now(),
            "patient": {"fasta_file_id": by_role["PATIENT_FASTA"], "metadata_file_id": by_role["PATIENT_METADATA"], "sequence_ids": list(sequence_ids)},
            "evidence": {"normalized_evidence_file_id": by_role["NORMALIZED_EVIDENCE"], "summary_file_id": by_role["EVIDENCE_SUMMARY"],
                         "source_registry_file_id": by_role["SOURCE_REGISTRY"], "validation_report_file_id": by_role["NORMALIZATION_VALIDATION_REPORT"]},
            "native_sources": native_sources, "alignments": [x.file_id for x in files if x.role == "ALIGNMENT_ARTIFACT"],
            "reads": [x.file_id for x in files if x.role.startswith("PREPARED_MINI_BAM")],
            "provenance": {"pipeline_manifest_file_id": by_role["PIPELINE_MANIFEST"], "tools_file_id": by_role["TOOL_PROVENANCE"],
                           "stages_file_id": by_role["STAGE_PROVENANCE"], "source_files_file_id": by_role["SOURCE_FILE_PROVENANCE"]},
            "config": {"run_config_file_id": by_role["RUN_CONFIG"], "locus_config_file_id": by_role["LOCUS_CONFIG"]},
            "files": [asdict(item) for item in files], "warnings": list(evidence.get("normalization_warnings", [])),
            "limitations": ["No consensus call or preferred caller is selected.", "No parental origin, expansion classification, or clinical interpretation is provided."]}
        atomic_write_json(work / "case-manifest.json", manifest)
        inventory_paths = sorted([item.relative_path for item in files] + ["case-manifest.json"])
        atomic_write_json(work / "checksums/sha256sums.json", {"record_schema_version": "1.0", "package_id": package_id,
            "excluded": ["checksums/sha256sums.json", "package-validation.json"],
            "files": [{"relative_path": rel, "sha256": sha256_file(work / rel), "size_bytes": (work / rel).stat().st_size} for rel in inventory_paths]})
        from .case_package_validation import validate_case_package
        validation = validate_case_package(work, write_report=True)
        if not validation["valid"]:
            details = "; ".join(f"{item.get('code')}: {item.get('message')}" for item in validation.get("issues", []))
            raise ConfigurationError(f"constructed case package failed validation: {details}")
        if destination.exists():
            backup = destination.with_name(f".{destination.name}.old-{uuid.uuid4().hex}")
            os.replace(destination, backup)
        os.replace(work, destination)
        if backup: shutil.rmtree(backup)
        return manifest
    except BaseException:
        if backup and backup.exists() and not destination.exists(): os.replace(backup, destination)
        shutil.rmtree(work, ignore_errors=True)
        raise


def validate_case_package(package_root: str | Path) -> dict[str, Any]:
    """Compatibility validator for the pre-Task-8 multi-locus manifest.

    New portable packages are validated by :mod:`case_package_validation`; the
    old public import remains stable for Task 1 contract consumers.
    """
    from .config import schema_path
    from .schema_validation import SchemaViolation, validate
    root = Path(package_root).resolve()
    manifest_path = root / "case-manifest.json" if root.is_dir() else root
    root = manifest_path.parent.resolve(); document = _json(manifest_path)
    if "package_contract_version" in document:
        from .case_package_validation import validate_case_package as independent
        return independent(root)
    try:
        validate(document, _json(schema_path("case-manifest.schema.json")))
    except SchemaViolation as exc:
        raise ConfigurationError(str(exc)) from exc
    files = document["files"]; file_ids = [item["file_id"] for item in files]
    if len(file_ids) != len(set(file_ids)): raise ConfigurationError("duplicate file_id in package inventory")
    known = set(file_ids); sequences: list[str] = []; references: list[str] = []
    for locus in document["loci"]:
        references.append(next((item["file_id"] for item in files if item["path"] == locus["locus_config"]["path"]), ""))
        for sequence in locus["patient_sequences"]:
            sequences.append(sequence["sequence_id"]); references.append(sequence["source_fasta_file_id"])
        references.extend(item["file_id"] for item in locus["native_caller_outputs"])
        references.extend(item["source_file_id"] for item in locus["normalized_evidence"])
    if len(sequences) != len(set(sequences)): raise ConfigurationError("duplicate sequence_id in case manifest")
    unknown = sorted(set(references) - known)
    if unknown: raise ConfigurationError("references to unknown file IDs: " + ", ".join(unknown))
    for item in files:
        posix = _safe_relative(item["path"]); candidate = root.joinpath(*posix.parts)
        if candidate.exists() and not candidate.resolve().is_relative_to(root): raise ConfigurationError(f"package path escapes through symlink: {item['path']}")
        if item["required"] and item["status"] == "AVAILABLE" and not candidate.is_file(): raise ConfigurationError(f"required available file is absent: {item['path']}")
    return document
