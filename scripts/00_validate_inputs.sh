#!/usr/bin/env bash
set -euo pipefail
# Validate configuration and immutable source inputs; callers are warning-only in scaffold mode.
usage(){ echo "Usage: $0 --config FILE [--dry-run]"; }
config=""; dry=false
while (($#)); do case "$1" in --config) config=${2:-}; shift 2;; --dry-run) dry=true; shift;; --help|-h) usage; exit 0;; *) echo "Error: unknown argument: $1" >&2; usage >&2; exit 2;; esac; done
[[ -n "$config" ]] || { echo "Error: --config is required" >&2; exit 2; }
[[ -f "$config" ]] || { echo "Error: configuration not found: $config" >&2; exit 1; }
python - "$config" "$dry" <<'PY'
import shutil, sys
from tr_calling_pipeline.config import load_config
cfg = load_config(sys.argv[1], check_inputs=True)
for name, executable in cfg.get("tools", {}).items():
    status = shutil.which(executable)
    print(f"Tool {name}: {status or 'WARNING: not found (allowed in scaffold development)'}")
print("Input validation passed" + (" (dry-run)" if sys.argv[2] == "true" else ""))
PY
