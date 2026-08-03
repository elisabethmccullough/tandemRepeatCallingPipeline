"""Small, dependency-free FASTA validation and deterministic writing."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import os
from pathlib import Path
import tempfile
from typing import Iterable


IUPAC_DNA = frozenset("ACGTRYSWKMBDHVN")


class FastaValidationError(ValueError):
    """Raised when a FASTA cannot be represented without guessing or data loss."""


@dataclass(frozen=True)
class FastaRecord:
    identifier: str
    description: str
    sequence: str


@dataclass(frozen=True)
class FastaRecordMetadata:
    record_schema_version: str
    record_id: str
    sequence_length: int
    sequence_sha256: str


def sequence_sha256(sequence: str) -> str:
    """Hash uppercase, unwrapped IUPAC letters (case is not biological identity)."""
    return hashlib.sha256(sequence.upper().encode("ascii")).hexdigest()


def read_fasta(path: str | Path) -> tuple[FastaRecord, ...]:
    source = Path(path)
    if not source.is_file() or source.stat().st_size == 0:
        raise FastaValidationError(f"FASTA is missing or empty: {source}")
    records: list[FastaRecord] = []
    identifier: str | None = None
    description = ""
    pieces: list[str] = []

    def finish() -> None:
        if identifier is None:
            return
        sequence = "".join(pieces)
        if not sequence:
            raise FastaValidationError(f"FASTA record {identifier!r} has an empty sequence")
        invalid = sorted(set(sequence.upper()) - IUPAC_DNA)
        if invalid:
            raise FastaValidationError(
                f"FASTA record {identifier!r} contains unsupported character(s): {''.join(invalid)}"
            )
        records.append(FastaRecord(identifier, description, sequence))

    try:
        lines = source.read_text(encoding="ascii").splitlines()
    except (UnicodeDecodeError, OSError) as exc:
        raise FastaValidationError(f"FASTA is not readable ASCII: {source}") from exc
    for line_number, line in enumerate(lines, 1):
        if line.startswith(">"):
            finish()
            header = line[1:].strip()
            if not header:
                raise FastaValidationError(f"empty FASTA identifier at line {line_number}")
            identifier, _, description = header.partition(" ")
            pieces = []
        else:
            if identifier is None:
                if not line.strip():
                    continue
                raise FastaValidationError(f"sequence before first FASTA header at line {line_number}")
            # Whitespace is formatting, but all non-whitespace symbols are biological content.
            pieces.append("".join(line.split()))
    finish()
    if not records:
        raise FastaValidationError(f"no FASTA records found: {source}")
    ids = [record.identifier for record in records]
    duplicates = sorted({item for item in ids if ids.count(item) > 1})
    if duplicates:
        raise FastaValidationError(f"duplicate FASTA record identifier(s): {', '.join(duplicates)}")
    return tuple(records)


def record_metadata(records: Iterable[FastaRecord]) -> tuple[FastaRecordMetadata, ...]:
    return tuple(FastaRecordMetadata("1.0", r.identifier, len(r.sequence), sequence_sha256(r.sequence)) for r in records)


def validation_report(path: str | Path) -> dict[str, object]:
    from .provenance import sha256_file

    source = Path(path).resolve()
    records = read_fasta(source)
    return {
        "schema_version": "1.0",
        "source_path": str(source),
        "source_file_sha256": sha256_file(source),
        "supported_alphabet": "ACGTRYSWKMBDHVN (case-insensitive; preserved as stored)",
        "sequence_hash_representation": "uppercase unwrapped IUPAC sequence bytes",
        "records": [asdict(item) for item in record_metadata(records)],
        "valid": True,
        "warnings": [],
    }


def write_fasta(path: str | Path, records: Iterable[FastaRecord], *, width: int = 80) -> None:
    """Atomically write records in supplied order, preserving sequence letter case."""
    if width < 1:
        raise ValueError("FASTA line width must be positive")
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    chunks: list[str] = []
    for record in records:
        if not record.identifier or any(c.isspace() for c in record.identifier):
            raise FastaValidationError(f"invalid output FASTA identifier: {record.identifier!r}")
        chunks.append(f">{record.identifier}\n")
        chunks.extend(record.sequence[i : i + width] + "\n" for i in range(0, len(record.sequence), width))
    fd, temporary = tempfile.mkstemp(prefix=f".{destination.name}.", dir=destination.parent)
    try:
        with os.fdopen(fd, "w", encoding="ascii", newline="\n") as handle:
            handle.writelines(chunks)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise
    read_fasta(destination)
