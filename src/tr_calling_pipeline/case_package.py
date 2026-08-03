"""Validation for portable GUI-facing case packages."""

from __future__ import annotations

import json
from pathlib import Path, PurePosixPath
from typing import Any

from .config import ConfigurationError, schema_path
from .schema_validation import SchemaViolation, validate


def _load_manifest(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ConfigurationError(f"Cannot read case manifest: {exc}") from exc
    if not isinstance(value, dict):
        raise ConfigurationError("case manifest root must be an object")
    return value


def validate_case_package(package_root: str | Path) -> dict[str, Any]:
    """Validate schema, referential integrity, and containment of inventory files."""
    root = Path(package_root).resolve()
    manifest_path = root / "case-manifest.json" if root.is_dir() else root
    root = manifest_path.parent.resolve()
    document = _load_manifest(manifest_path)
    schema_file = schema_path("case-manifest.schema.json")
    schema = json.loads(schema_file.read_text(encoding="utf-8"))
    try:
        validate(document, schema)
    except SchemaViolation as exc:
        raise ConfigurationError(str(exc)) from exc

    files = document["files"]
    file_ids = [item["file_id"] for item in files]
    if len(file_ids) != len(set(file_ids)):
        raise ConfigurationError("duplicate file_id in package inventory")
    known_files = set(file_ids)
    sequence_ids: list[str] = []
    references: list[str] = []
    for locus in document["loci"]:
        references.append(next((f["file_id"] for f in files if f["path"] == locus["locus_config"]["path"]), ""))
        for sequence in locus["patient_sequences"]:
            sequence_ids.append(sequence["sequence_id"])
            references.append(sequence["source_fasta_file_id"])
        references.extend(output["file_id"] for output in locus["native_caller_outputs"])
        references.extend(record["source_file_id"] for record in locus["normalized_evidence"])
    if len(sequence_ids) != len(set(sequence_ids)):
        raise ConfigurationError("duplicate sequence_id in case manifest")
    unknown = sorted(set(references) - known_files)
    if unknown:
        raise ConfigurationError(f"references to unknown file IDs: {', '.join(unknown)}")
    for item in files:
        posix = PurePosixPath(item["path"])
        if posix.is_absolute() or ".." in posix.parts or "\\" in item["path"]:
            raise ConfigurationError(f"unsafe package path: {item['path']}")
        candidate = root.joinpath(*posix.parts)
        if candidate.exists() and not candidate.resolve().is_relative_to(root):
            raise ConfigurationError(f"package path escapes through symlink: {item['path']}")
        if item["required"] and item["status"] == "AVAILABLE" and not candidate.is_file():
            raise ConfigurationError(f"required available file is absent: {item['path']}")
    return document
