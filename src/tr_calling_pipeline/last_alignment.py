"""Development-gated LAST preparation primitives.

No LAST release is claimed as laboratory verified.  Command spelling is kept in
one provisional adapter and is never selected without explicit opt-in.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from pathlib import Path
import re

from .provenance import sha256_file
from .tools import Tool, ToolId, detect_version, resolve_tool


class LastVersionClassification(str, Enum):
    VERIFIED_SUPPORTED = "VERIFIED_SUPPORTED"
    PROVISIONAL_DEVELOPMENT = "PROVISIONAL_DEVELOPMENT"
    UNSUPPORTED = "UNSUPPORTED"
    UNDETERMINED = "UNDETERMINED"


class UnsupportedLastVersion(ValueError): pass
class UnsupportedLastAlignment(ValueError): pass


@dataclass(frozen=True)
class LastCapabilities:
    alignment_format: str = "MAF"
    database_signature_suffixes: tuple[str, ...] = (".prj",)
    lastal_stdout_is_alignment: bool = True


@dataclass(frozen=True)
class LastCommandPlan:
    command_id: str
    argv: tuple[str, ...]
    stdout_path: str | None = None


@dataclass(frozen=True)
class LastAdapter:
    name: str = "last-development-provisional"
    capabilities: LastCapabilities = LastCapabilities()

    def lastdb_plan(self, executable: str, database_prefix: Path, fasta: Path,
                    additional_arguments=()) -> LastCommandPlan:
        return LastCommandPlan("lastdb", (executable, *tuple(additional_arguments), str(database_prefix), str(fasta)))

    def lastal_plan(self, executable: str, database_prefix: Path, reference: Path,
                    alignment: Path, additional_arguments=()) -> LastCommandPlan:
        return LastCommandPlan("lastal", (executable, *tuple(additional_arguments), str(database_prefix), str(reference)), str(alignment))


def classify_version(version: str | None) -> LastVersionClassification:
    if version is None: return LastVersionClassification.UNDETERMINED
    # LAST uses date-like numeric releases. Syntax remains explicitly provisional.
    if re.fullmatch(r"[0-9]{3,5}", version): return LastVersionClassification.PROVISIONAL_DEVELOPMENT
    return LastVersionClassification.UNSUPPORTED


def select_adapter(version: str | None, *, allow_provisional: bool = False) -> LastAdapter:
    kind = classify_version(version)
    if kind is LastVersionClassification.PROVISIONAL_DEVELOPMENT and allow_provisional: return LastAdapter()
    if kind is LastVersionClassification.UNDETERMINED:
        raise UnsupportedLastVersion("LAST version could not be determined; syntax will not be guessed")
    if kind is LastVersionClassification.PROVISIONAL_DEVELOPMENT:
        raise UnsupportedLastVersion("LAST has no laboratory-verified adapter; provisional adapter is disabled")
    raise UnsupportedLastVersion(f"unsupported LAST version: {version}")


def detect_last(executable: str, config_directory: Path, required: bool, tool_id: ToolId) -> Tool:
    if tool_id not in (ToolId.LASTDB, ToolId.LASTAL): raise ValueError("tool_id must be LASTDB or LASTAL")
    tool = resolve_tool(Tool(tool_id, tool_id.value.lower(), executable, required), config_directory)
    return detect_version(tool, pattern=r"(?i)(?:last(?:db|al)?[^0-9]*)?([0-9]{3,5})") if tool.resolved_executable else tool


def write_single_record_fasta(path: str | Path, sequence_id: str, sequence: str, width: int = 80) -> str:
    """Write deterministic FASTA without altering sequence letters."""
    if not sequence or "\n" in sequence or "\r" in sequence: raise ValueError("sequence must be non-empty and newline-free")
    target = Path(path); target.parent.mkdir(parents=True, exist_ok=True)
    text = f">{sequence_id}\n" + "\n".join(sequence[i:i+width] for i in range(0, len(sequence), width)) + "\n"
    target.write_text(text, encoding="ascii", newline="\n")
    return sha256_file(target)


def validate_maf(path: str | Path, expected_query_id: str) -> tuple[str, ...]:
    source = Path(path)
    if not source.is_file(): raise UnsupportedLastAlignment("alignment is missing")
    if source.stat().st_size == 0: raise UnsupportedLastAlignment("alignment is empty")
    lines = source.read_text(encoding="utf-8").splitlines()
    if not lines or not lines[0].startswith("##maf"): raise UnsupportedLastAlignment("unsupported alignment format")
    blocks = [i for i, line in enumerate(lines) if line.startswith("a ") or line == "a"]
    if not blocks: raise UnsupportedLastAlignment("MAF has no alignment block")
    identifiers = tuple(line.split()[1] for line in lines if line.startswith("s ") and len(line.split()) >= 7)
    if expected_query_id not in identifiers: raise UnsupportedLastAlignment("expected patient sequence is absent from MAF")
    return identifiers
