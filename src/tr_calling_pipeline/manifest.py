"""Run-manifest creation and incremental persistence."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import subprocess
from typing import Any, Iterable


def file_checksum(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def repository_commit(cwd: str | Path | None = None) -> str | None:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=cwd, text=True, capture_output=True, check=False
    )
    return result.stdout.strip() if result.returncode == 0 else None


def create_manifest(config: dict[str, Any], config_path: str | Path) -> dict[str, Any]:
    now = datetime.now(timezone.utc)
    inputs = dict(config["inputs"])
    resolved_paths = dict(inputs)
    if "locus_config" in config:
        resolved_paths["locus_config"] = config["locus_config"]
    if "output_root" in config["run"]:
        resolved_paths["output_root"] = config["run"]["output_root"]
    checksums = {
        key: file_checksum(path) for key, path in inputs.items() if Path(path).is_file()
    }
    return {
        "run_id": f'{config["run"]["sample_id"]}_{config["run"]["locus_id"]}_{now.strftime("%Y%m%dT%H%M%SZ")}',
        "sample_id": config["run"]["sample_id"], "locus_id": config["run"]["locus_id"],
        "start_time": now.isoformat(), "repository_commit": repository_commit(),
        "configuration_path": str(Path(config_path).resolve()), "resolved_inputs": inputs,
        "resolved_paths": resolved_paths,
        "input_checksums": checksums, "tool_versions": {}, "stages": {},
        "output_paths": [], "warnings": [], "completion_status": "running",
    }


def write_manifest(manifest: dict[str, Any], path: str | Path) -> Path:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return destination


def update_stage(manifest: dict[str, Any], stage: str, status: str, outputs: Iterable[str] = ()) -> None:
    manifest["stages"][stage] = {"status": status, "updated_at": datetime.now(timezone.utc).isoformat()}
    manifest["output_paths"].extend(str(item) for item in outputs)
