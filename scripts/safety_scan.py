#!/usr/bin/env python3
"""Repository checks that are deterministic and dependency-free."""
from pathlib import Path
import re

root = Path(__file__).resolve().parents[1]
errors = []
for path in root.rglob("*"):
    if ".git" in path.parts or not path.is_file():
        continue
    relative = path.relative_to(root)
    if path.is_symlink(): errors.append(f"tracked/scanned symlink: {relative}")
    if path.stat().st_size > 5_000_000: errors.append(f"large file: {relative}")
    if path.suffix == ".py":
        text = path.read_text(encoding="utf-8", errors="replace")
        if re.search(r"\beval\s*\(", text): errors.append(f"eval use: {relative}")
        unsafe_shell_pattern = r"shell\s*=\s*" + "True"
        if re.search(unsafe_shell_pattern, text): errors.append(f"unsafe shell option: {relative}")
for error in errors: print(error)
raise SystemExit(bool(errors))
