"""Verification vocabulary and repository/resource contract checks."""

from __future__ import annotations

import json
from importlib.resources import files
from pathlib import Path
from typing import Iterable

VERIFICATION_LEVELS = (
    "UNIT_TESTED", "SYNTHETIC_INTEGRATION_TESTED", "REAL_TOOL_SMOKE_TESTED",
    "LABORATORY_VERIFIED", "UNVERIFIED",
)


def schema_directory() -> Path:
    """Return installed schemas, falling back to a source checkout for developers."""
    packaged = files("tr_calling_pipeline").joinpath("schemas")
    if packaged.is_dir():
        return Path(str(packaged))
    checkout = Path(__file__).resolve().parents[2] / "schemas"
    if checkout.is_dir():
        return checkout
    raise FileNotFoundError("no packaged schema directory is available")


def _refs(value: object) -> Iterable[str]:
    if isinstance(value, dict):
        for key, child in value.items():
            if key == "$ref" and isinstance(child, str):
                yield child
            yield from _refs(child)
    elif isinstance(value, list):
        for child in value:
            yield from _refs(child)


def validate_schemas(directory: Path | None = None) -> dict[str, object]:
    root = directory or schema_directory()
    errors: list[str] = []
    ids: dict[str, str] = {}
    paths = sorted(root.glob("*.schema.json"))
    for path in paths:
        try:
            document = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            errors.append(f"{path.name}: {exc}")
            continue
        schema_id = document.get("$id")
        if schema_id:
            if schema_id in ids:
                errors.append(f"duplicate $id {schema_id!r}: {ids[schema_id]} and {path.name}")
            ids[schema_id] = path.name
        for ref in _refs(document):
            if ref.startswith(("http://", "https://")):
                errors.append(f"{path.name}: network $ref is prohibited: {ref}")
            elif not ref.startswith("#") and not (root / ref.split("#", 1)[0]).is_file():
                errors.append(f"{path.name}: missing local $ref target: {ref}")
    return {"valid": not errors, "schema_count": len(paths), "errors": errors}
