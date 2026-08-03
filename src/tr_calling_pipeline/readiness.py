"""Conservative release-readiness reporting."""
from __future__ import annotations
import os
from pathlib import Path
from .verification import validate_schemas

def check_release_readiness(root: Path | None = None) -> tuple[str, list[str]]:
    root = root or Path.cwd()
    failures = []
    if not validate_schemas()["valid"]: failures.append("schema validation failed")
    for name in ("README.md", "CHANGELOG.md", "LICENSE"):
        if not (root / name).is_file(): failures.append(f"missing {name}")
    evidence = {
        "full synthetic runner demo": os.environ.get("TR_PIPELINE_FULL_DEMO_VERIFIED") == "1",
        "clean wheel installation": os.environ.get("TR_PIPELINE_WHEEL_VERIFIED") == "1",
        "hosted CI matrix": os.environ.get("TR_PIPELINE_HOSTED_CI_VERIFIED") == "1",
    }
    failures.extend(f"missing successful evidence: {name}" for name, present in evidence.items() if not present)
    if failures: return "NOT_READY", failures
    return "READY_WITH_LIMITATIONS", ["real caller installations and laboratory workflows are not verified"]
