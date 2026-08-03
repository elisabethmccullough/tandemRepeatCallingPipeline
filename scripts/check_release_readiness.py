#!/usr/bin/env python3
"""Run the lightweight, non-claiming release-readiness report."""
from pathlib import Path
from tr_calling_pipeline.readiness import check_release_readiness

status, details = check_release_readiness(Path.cwd())
print(status)
for detail in details:
    print(f"- {detail}")
raise SystemExit(2 if status == "NOT_READY" else 0)
