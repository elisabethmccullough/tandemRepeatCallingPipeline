"""Command-line interface for shared contract validation."""

from __future__ import annotations

import argparse
import json
from typing import Sequence

from .case_package import validate_case_package
from .config import load_config, load_locus_config


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="tr-pipeline")
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("validate-run-config", "validate-locus-config", "validate-case-package", "print-resolved-config"):
        command = commands.add_parser(name)
        command.add_argument("path")
    args = parser.parse_args(argv)
    if args.command == "validate-run-config":
        load_config(args.path)
    elif args.command == "validate-locus-config":
        load_locus_config(args.path)
    elif args.command == "validate-case-package":
        validate_case_package(args.path)
    else:
        print(json.dumps(load_config(args.path), indent=2, sort_keys=True))
    return 0
