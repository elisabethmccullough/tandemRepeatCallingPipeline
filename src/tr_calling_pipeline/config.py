"""Configuration loading, validation, and path resolution."""

from __future__ import annotations

from pathlib import Path
import re
from typing import Any

import yaml

SAFE_IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
REQUIRED_FIELDS = (
    "run.sample_id", "run.locus_id", "run.output_root",
    "inputs.assembly_fasta", "inputs.mini_bam", "inputs.mini_bam_index",
    "inputs.reference_fasta", "inputs.reference_fasta_index", "locus_config",
    "execution.threads", "execution.overwrite",
)


class ConfigurationError(ValueError):
    """Raised when pipeline configuration is incomplete or unsafe."""


def _get(data: dict[str, Any], dotted: str) -> Any:
    value: Any = data
    for part in dotted.split("."):
        if not isinstance(value, dict) or part not in value:
            raise ConfigurationError(f"Missing required configuration field: {dotted}")
        value = value[part]
    return value


def validate_identifier(value: str, field: str = "identifier") -> str:
    if not isinstance(value, str) or not SAFE_IDENTIFIER.fullmatch(value):
        raise ConfigurationError(
            f"{field} must contain only letters, numbers, '.', '_' or '-' and cannot be empty"
        )
    return value


def load_config(path: str | Path, *, check_inputs: bool = False) -> dict[str, Any]:
    """Load YAML, validate its model, and resolve paths from the repository root.

    Relative paths are interpreted from the current working directory, allowing the
    same checked-in configuration to be invoked consistently from the repository.
    """
    config_path = Path(path).expanduser()
    if not config_path.is_file():
        raise ConfigurationError(f"Configuration file does not exist: {config_path}")
    try:
        data = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise ConfigurationError(f"Invalid YAML in {config_path}: {exc}") from exc
    if not isinstance(data, dict):
        raise ConfigurationError("Configuration root must be a YAML mapping")
    for field in REQUIRED_FIELDS:
        _get(data, field)
    validate_identifier(_get(data, "run.sample_id"), "run.sample_id")
    validate_identifier(_get(data, "run.locus_id"), "run.locus_id")
    if not isinstance(_get(data, "execution.threads"), int) or data["execution"]["threads"] < 1:
        raise ConfigurationError("execution.threads must be a positive integer")

    base = Path.cwd()
    data["_config_path"] = str(config_path.resolve())
    for key, value in data["inputs"].items():
        resolved = (base / value).resolve() if not Path(value).is_absolute() else Path(value).resolve()
        data["inputs"][key] = str(resolved)
        if check_inputs and (not resolved.is_file() or resolved.stat().st_size == 0):
            raise ConfigurationError(f"Input file is missing or empty ({key}): {resolved}")
    for field in ("locus_config",):
        value = Path(data[field]).expanduser()
        data[field] = str(((base / value) if not value.is_absolute() else value).resolve())
        if check_inputs and not Path(data[field]).is_file():
            raise ConfigurationError(f"Locus configuration does not exist: {data[field]}")
    output = Path(data["run"]["output_root"]).expanduser()
    data["run"]["output_root"] = str(((base / output) if not output.is_absolute() else output).resolve())
    return data
