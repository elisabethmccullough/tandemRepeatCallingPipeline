"""Schema-backed configuration loading, semantic validation, and path resolution."""

from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import re
from typing import Any

import yaml
from .schema_validation import SchemaViolation, validate

SAFE_IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
PATH_INPUTS = ("assembly_fasta", "mini_bam", "mini_bam_index", "reference_fasta", "reference_fasta_index")


class ConfigurationError(ValueError):
    """Raised when pipeline configuration is incomplete or unsafe."""


def schema_path(name: str) -> Path:
    return Path(__file__).resolve().parents[2] / "schemas" / name


def _validate_schema(data: Any, name: str) -> None:
    schema = json.loads(schema_path(name).read_text(encoding="utf-8"))
    try:
        validate(data, schema)
    except SchemaViolation as exc:
        raise ConfigurationError(str(exc)) from exc


def validate_identifier(value: str, field: str = "identifier") -> str:
    if not isinstance(value, str) or not SAFE_IDENTIFIER.fullmatch(value):
        raise ConfigurationError(f"{field} must contain only letters, numbers, '.', '_' or '-' and cannot be empty")
    return value


def _load_yaml(path: str | Path) -> tuple[dict[str, Any], Path]:
    source = Path(path).expanduser().resolve()
    if not source.is_file():
        raise ConfigurationError(f"Configuration file does not exist: {source}")
    try:
        data = yaml.safe_load(source.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise ConfigurationError(f"Invalid YAML in {source}: {exc}") from exc
    if not isinstance(data, dict):
        raise ConfigurationError("Configuration root must be a YAML mapping")
    return data, source


def load_config(path: str | Path, *, check_inputs: bool = False) -> dict[str, Any]:
    """Validate a run YAML and resolve all paths from that YAML's directory."""
    original, config_path = _load_yaml(path)
    _validate_schema(original, "run-config.schema.json")
    data = deepcopy(original)
    records = data["inputs"]["assembly_records"]
    sequence_ids = [record["sequence_id"] for record in records]
    record_ids = [record["record_id"] for record in records]
    if len(sequence_ids) != len(set(sequence_ids)) or len(record_ids) != len(set(record_ids)):
        raise ConfigurationError("assembly record_id and sequence_id values must be unique")
    if any(r["sequence_role"] not in {"PATIENT_HAPLOTYPE", "PATIENT_CONSENSUS"} for r in records):
        raise ConfigurationError("patient assembly records may only use PATIENT_HAPLOTYPE or PATIENT_CONSENSUS")

    base = config_path.parent
    data["_config_path"] = str(config_path)
    for key in PATH_INPUTS:
        candidate = Path(data["inputs"][key]).expanduser()
        resolved = candidate if candidate.is_absolute() else (base / candidate).resolve()
        data["inputs"][key] = str(resolved)
        if check_inputs and (not resolved.is_file() or resolved.stat().st_size == 0):
            raise ConfigurationError(f"Input file is missing or empty ({key}): {resolved}")
    for container, key in ((data, "locus_config"), (data["run"], "output_root"), (data["case_package"], "package_root")):
        candidate = Path(container[key]).expanduser()
        container[key] = str(candidate if candidate.is_absolute() else (base / candidate).resolve())
    if check_inputs and not Path(data["locus_config"]).is_file():
        raise ConfigurationError(f"Locus configuration does not exist: {data['locus_config']}")
    return data


def load_locus_config(path: str | Path) -> dict[str, Any]:
    data, _ = _load_yaml(path)
    _validate_schema(data, "locus-config.schema.json")
    block_ids = [b["block_id"] for b in data["repeat_blocks"]]
    orders = [b["order"] for b in data["repeat_blocks"]]
    if len(block_ids) != len(set(block_ids)) or len(orders) != len(set(orders)):
        raise ConfigurationError("repeat block IDs and order values must be unique")
    motif_ids = [m["motif_id"] for b in data["repeat_blocks"] for m in b["motifs"]]
    if len(motif_ids) != len(set(motif_ids)):
        raise ConfigurationError("motif IDs must be unique across the locus")
    region = data["locus"]["target_region"]
    if data["release_status"] != "development" and (region["start"] is None or region["end"] is None):
        raise ConfigurationError("only development locus configurations may use null target coordinates")
    if region["start"] is not None and region["end"] is not None and region["start"] >= region["end"]:
        raise ConfigurationError("target region start must precede end")
    return data
