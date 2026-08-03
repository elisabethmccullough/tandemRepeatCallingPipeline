"""Independent validation of a moved case package using package contents only."""

from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path, PurePosixPath
from typing import Any

from .case_package import (CONTROL_FILES, CasePackageValidationIssue,
                           CasePackageValidationReport, _fasta, _json, _safe_relative)
from .provenance import atomic_write_json, sha256_file, utc_now
from .schema_validation import SchemaViolation, validate
from .config import schema_path


PACKAGED_JSON_SCHEMAS = {
    "case-manifest.json": "case-package-manifest.schema.json",
    "evidence/normalized-evidence.json": "unified-normalized-evidence.schema.json",
    "evidence/evidence-summary.json": "unified-evidence-summary.schema.json",
    "evidence/source-registry.json": "unified-source-registry.schema.json",
    "evidence/normalization-validation-report.json": "unified-validation-report.schema.json",
    "checksums/sha256sums.json": "case-package-checksums.schema.json",
    "package-validation.json": "case-package-validation.schema.json",
}


def validate_case_package(package_root: str | Path, *, write_report: bool = False) -> dict[str, Any]:
    root = Path(package_root).resolve()
    if root.is_file(): root = root.parent
    issues: list[CasePackageValidationIssue] = []
    validated: list[str] = []

    def issue(severity: str, code: str, message: str, relative_path: str | None = None,
              file_id: str | None = None, record_id: str | None = None,
              caller: str | None = None, associated_sequence_id: str | None = None) -> None:
        issues.append(CasePackageValidationIssue(f"issue-{len(issues)+1:04d}", severity, code, message,
            relative_path, file_id, record_id, caller, associated_sequence_id))

    documents: dict[str, dict[str, Any]] = {}
    for relative, schema_name in PACKAGED_JSON_SCHEMAS.items():
        path = root / relative
        # The builder invokes validation before creating the first validation
        # report. All independently validated, published packages require it.
        if relative == "package-validation.json" and write_report and not path.exists():
            continue
        try:
            document = _json(path)
            documents[relative] = document
            schema = json.loads(schema_path(schema_name).read_text(encoding="utf-8"))
            validate(document, schema)
        except Exception as exc:
            issue("ERROR", "SCHEMA_VALIDATION_ERROR", f"{schema_name}: {exc}", relative)
    manifest = documents.get("case-manifest.json", {})
    files = manifest.get("files", []) if isinstance(manifest.get("files", []), list) else []
    ids: set[str] = set(); paths: set[str] = set()
    file_index: dict[str, dict[str, Any]] = {}
    for entry in files:
        file_id, relative = entry.get("file_id"), entry.get("relative_path")
        if file_id in ids: issue("ERROR", "DUPLICATE_FILE_ID", f"duplicate file ID: {file_id}", relative, file_id)
        if relative in paths: issue("ERROR", "DUPLICATE_RELATIVE_PATH", f"duplicate path: {relative}", relative, file_id)
        ids.add(file_id); paths.add(relative)
        if file_id not in file_index: file_index[file_id] = entry
        try: posix = _safe_relative(relative)
        except Exception as exc:
            issue("ERROR", "MANIFEST_PATH_TRAVERSAL", str(exc), relative, file_id); continue
        candidate = root.joinpath(*posix.parts)
        if candidate.is_symlink():
            issue("ERROR", "UNSAFE_SYMLINK", "registered package paths may not be symlinks", relative, file_id); continue
        if not candidate.is_file():
            issue("ERROR", "MISSING_REQUIRED_FILE" if entry.get("required") else "MISSING_FILE", "registered file is absent", relative, file_id); continue
        if not candidate.resolve().is_relative_to(root):
            issue("ERROR", "PATH_ESCAPE", "registered file escapes package", relative, file_id); continue
        actual_size, actual_hash = candidate.stat().st_size, sha256_file(candidate)
        if actual_size != entry.get("size_bytes"): issue("ERROR", "SIZE_MISMATCH", "file size differs from manifest", relative, file_id)
        if actual_hash != entry.get("sha256"): issue("ERROR", "CHECKSUM_MISMATCH", "file checksum differs from manifest", relative, file_id)
        validated.append(relative)

    checksum_paths: set[str] = set()
    try:
        checksum = documents.get("checksums/sha256sums.json") or _json(root / "checksums/sha256sums.json")
        for entry in checksum.get("files", []):
            relative = entry["relative_path"]
            if relative in checksum_paths: issue("ERROR", "DUPLICATE_CHECKSUM_PATH", "duplicate checksum inventory path", relative)
            checksum_paths.add(relative)
            posix = _safe_relative(relative); candidate = root.joinpath(*posix.parts)
            if not candidate.is_file(): issue("ERROR", "CHECKSUM_FILE_MISSING", "checksum entry is absent", relative); continue
            if candidate.stat().st_size != entry["size_bytes"]: issue("ERROR", "CHECKSUM_SIZE_MISMATCH", "checksum inventory size differs", relative)
            if sha256_file(candidate) != entry["sha256"]: issue("ERROR", "CHECKSUM_MISMATCH", "checksum inventory digest differs", relative)
    except Exception as exc:
        issue("ERROR", "INVALID_CHECKSUM_INVENTORY", str(exc), "checksums/sha256sums.json")
    expected_checksums = paths | {"case-manifest.json"}
    for relative in sorted(expected_checksums - checksum_paths): issue("ERROR", "MISSING_CHECKSUM_ENTRY", "file absent from checksum inventory", relative)
    for relative in sorted(checksum_paths - expected_checksums): issue("ERROR", "UNREGISTERED_CHECKSUM_ENTRY", "unexpected checksum entry", relative)

    actual: set[str] = set()
    for candidate in root.rglob("*") if root.is_dir() else ():
        relative = candidate.relative_to(root).as_posix()
        if candidate.is_symlink(): issue("ERROR", "UNSAFE_SYMLINK", "package contains a symlink", relative)
        elif candidate.is_file(): actual.add(relative)
    for relative in sorted(actual - paths - CONTROL_FILES): issue("ERROR", "UNEXPECTED_UNREGISTERED_FILE", "file is not registered", relative)

    sequence_ids = set(manifest.get("patient", {}).get("sequence_ids", []))
    metadata_sequence_ids: set[str] = set()
    try:
        fasta = _fasta(root / "patient/patient-sequences.fasta")
        metadata = _json(root / "patient/patient-sequences.metadata.json")
        metadata_sequence_ids = {item["sequence_id"] for item in metadata["sequences"]}
        if len(metadata_sequence_ids) != len(metadata["sequences"]): issue("ERROR", "DUPLICATE_PATIENT_SEQUENCE_ID", "patient metadata contains duplicate sequence IDs")
        if set(fasta) != metadata_sequence_ids or metadata_sequence_ids != sequence_ids:
            issue("ERROR", "PATIENT_SEQUENCE_SET_MISMATCH", "FASTA, metadata, and manifest sequence IDs differ")
        import hashlib
        for item in metadata["sequences"]:
            if item["sequence_id"] in fasta and hashlib.sha256(fasta[item["sequence_id"]].encode("ascii")).hexdigest() != item["sequence_sha256"]:
                issue("ERROR", "PATIENT_FASTA_CHECKSUM_MISMATCH", "biological sequence checksum differs", associated_sequence_id=item["sequence_id"])
            if not item.get("source_fasta_record_id"): issue("ERROR", "MISSING_SOURCE_RECORD_ID", "source record ID is absent", associated_sequence_id=item["sequence_id"])
    except Exception as exc: issue("ERROR", "INVALID_PATIENT_DATA", str(exc))
    try:
        evidence = documents.get("evidence/normalized-evidence.json") or _json(root / "evidence/normalized-evidence.json")
        registry = documents.get("evidence/source-registry.json") or _json(root / "evidence/source-registry.json")
        source_counts: dict[str, int] = {}
        for source in registry.get("sources", []): source_counts[source["file_id"]] = source_counts.get(source["file_id"], 0) + 1
        source_ids = set(source_counts)
        for source_id, count in source_counts.items():
            if count != 1: issue("ERROR", "DUPLICATE_SOURCE_REGISTRY_FILE_ID", "source registry file ID does not resolve exactly once", file_id=source_id)
        registered_native = {item.get("source_file_id"): item for item in files if item.get("source_file_id")}
        for source in registry.get("sources", []):
            packaged = registered_native.get(source.get("file_id"))
            if packaged is None:
                # Native output inclusion is configurable; omission is not corruption.
                continue
            if packaged.get("sha256") != source.get("sha256") or packaged.get("size_bytes") != source.get("size_bytes"):
                issue("ERROR", "NATIVE_SOURCE_CHECKSUM_MISMATCH", "packaged native identity differs from source registry",
                      packaged.get("relative_path"), packaged.get("file_id"), caller=source.get("caller"),
                      associated_sequence_id=source.get("associated_sequence_id"))
        identity = ("case_id", "subject_id", "sample_id", "locus_id")
        for key in identity:
            if evidence.get(key) != manifest.get(key): issue("ERROR", "EVIDENCE_IDENTITY_MISMATCH", f"evidence {key} differs")
        for record in evidence.get("records", []):
            sequence = record.get("associated_sequence_id")
            if sequence is not None and sequence not in sequence_ids: issue("ERROR", "UNKNOWN_SEQUENCE_ASSOCIATION", "evidence refers to unknown sequence", record_id=record.get("record_id"), caller=record.get("caller"), associated_sequence_id=sequence)
            if source_counts.get(record.get("source_file_id"), 0) != 1: issue("ERROR", "UNKNOWN_SOURCE_FILE", "evidence source file does not resolve exactly once", record_id=record.get("record_id"))
        report = documents.get("evidence/normalization-validation-report.json") or _json(root / "evidence/normalization-validation-report.json")
        if report.get("valid") is not True: issue("ERROR", "NORMALIZATION_REPORT_INVALID", "packaged Stage 08 report is invalid")
    except Exception as exc: issue("ERROR", "INVALID_EVIDENCE", str(exc))

    def reference(file_id: Any, expected_roles: set[str], location: str) -> None:
        entry = file_index.get(file_id)
        if entry is None:
            issue("ERROR", "UNKNOWN_MANIFEST_FILE_ID", f"{location} references unknown file ID: {file_id}", file_id=file_id)
        elif entry.get("role") not in expected_roles:
            issue("ERROR", "WRONG_MANIFEST_FILE_ROLE", f"{location} requires {sorted(expected_roles)}, got {entry.get('role')}",
                  entry.get("relative_path"), file_id)

    patient = manifest.get("patient", {})
    reference(patient.get("fasta_file_id"), {"PATIENT_FASTA"}, "patient.fasta_file_id")
    reference(patient.get("metadata_file_id"), {"PATIENT_METADATA"}, "patient.metadata_file_id")
    evidence_refs = manifest.get("evidence", {})
    for field, role in (("normalized_evidence_file_id", "NORMALIZED_EVIDENCE"), ("summary_file_id", "EVIDENCE_SUMMARY"),
                        ("source_registry_file_id", "SOURCE_REGISTRY"), ("validation_report_file_id", "NORMALIZATION_VALIDATION_REPORT")):
        reference(evidence_refs.get(field), {role}, f"evidence.{field}")
    provenance_roles = (("pipeline_manifest_file_id", "PIPELINE_MANIFEST"), ("tools_file_id", "TOOL_PROVENANCE"),
                        ("stages_file_id", "STAGE_PROVENANCE"), ("source_files_file_id", "SOURCE_FILE_PROVENANCE"))
    for field, role in provenance_roles: reference(manifest.get("provenance", {}).get(field), {role}, f"provenance.{field}")
    for field, role in (("run_config_file_id", "RUN_CONFIG"), ("locus_config_file_id", "LOCUS_CONFIG")):
        reference(manifest.get("config", {}).get(field), {role}, f"config.{field}")
    for file_id in manifest.get("alignments", []): reference(file_id, {"ALIGNMENT_ARTIFACT"}, "alignments")
    for file_id in manifest.get("reads", []): reference(file_id, {"PREPARED_MINI_BAM", "PREPARED_MINI_BAM_INDEX"}, "reads")
    native_seen: set[Any] = set()
    registry_counts: dict[Any, int] = {}
    for item in documents.get("evidence/source-registry.json", {}).get("sources", []):
        registry_counts[item.get("file_id")] = registry_counts.get(item.get("file_id"), 0) + 1
    for native in manifest.get("native_sources", []):
        native_id = native.get("file_id")
        if native_id in native_seen: issue("ERROR", "DUPLICATE_NATIVE_SOURCE_REFERENCE", "duplicate native source reference", file_id=native_id)
        native_seen.add(native_id)
        reference(native.get("file_id"), {"CALLER_NATIVE_OUTPUT"}, "native_sources.file_id")
        if registry_counts.get(native.get("source_file_id"), 0) != 1:
            issue("ERROR", "UNKNOWN_NATIVE_SOURCE_FILE_ID", "native source_file_id does not resolve exactly once in packaged source registry",
                  file_id=native.get("file_id"))
    for entry in files:
        sequence = entry.get("associated_sequence_id")
        if sequence is not None and sequence not in metadata_sequence_ids:
            issue("ERROR", "UNKNOWN_FILE_SEQUENCE_ASSOCIATION", "manifest file refers to unknown patient sequence",
                  entry.get("relative_path"), entry.get("file_id"), associated_sequence_id=sequence)
    for relative in ("provenance/pipeline-run-manifest.json", "provenance/tools.json", "provenance/stages.json",
                     "provenance/source-files.json", "config/run-config.yaml", "config/locus-config.yaml"):
        if not (root / relative).is_file(): issue("ERROR", "MISSING_PROVENANCE", "required provenance/config file absent", relative)

    errors = sum(item.severity == "ERROR" for item in issues); warnings = sum(item.severity == "WARNING" for item in issues)
    report = CasePackageValidationReport("1.0", manifest.get("package_id"), not errors, utc_now(), errors, warnings,
                                         tuple(issues), tuple(sorted(set(validated))))
    value = report.to_dict()
    if write_report: atomic_write_json(root / "package-validation.json", value)
    return value
