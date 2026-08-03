"""Caller-neutral immutable registration of native evidence files."""
from __future__ import annotations
from dataclasses import asdict, dataclass
from pathlib import Path
import mimetypes
from .provenance import sha256_file

@dataclass(frozen=True)
class NativeCallerOutput:
    file_id: str
    role: str
    path: str
    media_type: str
    sha256: str
    size_bytes: int
    caller: str
    caller_version: str | None
    analysis_source: str
    producer_command_id: str

    @classmethod
    def from_path(cls, path: str | Path, *, file_id: str, caller_version: str | None,
                  analysis_source: str, producer_command_id: str) -> "NativeCallerOutput":
        source = Path(path)
        return cls(file_id, "NATIVE_CALLER_OUTPUT", str(source),
                   mimetypes.guess_type(source.name)[0] or "application/octet-stream",
                   sha256_file(source), source.stat().st_size, "VAMOS", caller_version,
                   analysis_source, producer_command_id)

    def to_dict(self):
        return asdict(self)
