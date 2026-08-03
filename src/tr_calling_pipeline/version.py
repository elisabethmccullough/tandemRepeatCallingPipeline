"""Pipeline package and source checkout identity."""
from dataclasses import dataclass
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
import subprocess

@dataclass(frozen=True)
class PipelineIdentity:
    pipeline_name: str; pipeline_version: str; git_commit: str | None; git_dirty: bool | None

def pipeline_identity(directory: str | Path | None = None) -> PipelineIdentity:
    try: package_version = version("tandem-repeat-calling-pipeline")
    except PackageNotFoundError: package_version = "0.1.0.dev0"
    cwd = Path(directory).resolve() if directory else Path(__file__).resolve().parents[2]
    try:
        commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=cwd, text=True, capture_output=True, check=True).stdout.strip()
        dirty = bool(subprocess.run(["git", "status", "--porcelain"], cwd=cwd, text=True, capture_output=True, check=True).stdout)
    except (OSError, subprocess.SubprocessError): commit, dirty = None, None
    return PipelineIdentity("tandem-repeat-calling-pipeline", package_version, commit, dirty)
