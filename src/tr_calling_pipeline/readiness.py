"""Conservative release-readiness reporting."""
from __future__ import annotations
from pathlib import Path
from .verification import validate_schemas

def check_release_readiness(root: Path | None = None) -> tuple[str, list[str]]:
    root = root or Path.cwd()
    failures = []
    if not validate_schemas()["valid"]: failures.append("schema validation failed")
    for name in ("README.md", "CHANGELOG.md", "LICENSE"):
        if not (root / name).is_file(): failures.append(f"missing {name}")
    if failures: return "NOT_READY", failures
    return "READY_WITH_LIMITATIONS", ["real caller installations and laboratory workflows are not verified"]
