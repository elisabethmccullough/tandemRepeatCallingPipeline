"""Deterministic and crash-safe provenance primitives."""
from __future__ import annotations
from dataclasses import asdict, dataclass, is_dataclass
from datetime import datetime, timezone
import hashlib, json, os, tempfile
from pathlib import Path
from typing import Any

def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="microseconds").replace("+00:00", "Z")

def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""): digest.update(block)
    return digest.hexdigest()

@dataclass(frozen=True)
class FileIdentity:
    path: str
    sha256: str
    size_bytes: int

def file_identity(path: str | Path) -> FileIdentity:
    p = Path(path).resolve()
    return FileIdentity(str(p), sha256_file(p), p.stat().st_size)

def canonical_digest(value: Any) -> str:
    def clean(v: Any) -> Any:
        if is_dataclass(v): v = asdict(v)
        if isinstance(v, dict): return {k: clean(x) for k, x in sorted(v.items()) if not k.startswith("_")}
        if isinstance(v, (list, tuple)): return [clean(x) for x in v]
        if isinstance(v, Path): return str(v)
        return v
    raw = json.dumps(clean(value), sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(raw.encode()).hexdigest()

def atomic_write_json(path: str | Path, value: Any) -> None:
    target = Path(path); target.parent.mkdir(parents=True, exist_ok=True)
    payload = asdict(value) if is_dataclass(value) else value
    fd, temporary = tempfile.mkstemp(prefix=f".{target.name}.", dir=target.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True); handle.write("\n"); handle.flush(); os.fsync(handle.fileno())
        os.replace(temporary, target)
    except BaseException:
        try: os.unlink(temporary)
        except FileNotFoundError: pass
        raise
