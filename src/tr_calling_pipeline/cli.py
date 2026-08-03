"""Command-line interface for shared contract validation."""

from __future__ import annotations

import argparse
import json
from typing import Sequence

from .case_package import validate_case_package
from .config import load_config, load_locus_config
from .runner import run


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="tr-pipeline")
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("validate-run-config", "validate-locus-config", "validate-case-package", "print-resolved-config"):
        command = commands.add_parser(name)
        command.add_argument("path")
    command = commands.add_parser("run")
    command.add_argument("--config", required=True)
    command.add_argument("--dry-run", action="store_true")
    command.add_argument("--start-stage")
    command.add_argument("--stop-stage")
    command.add_argument("--resume", action=argparse.BooleanOptionalAction, default=True)
    command.add_argument("--overwrite", action="store_true")
    command.add_argument("--execution-mode", choices=("NATIVE", "APPTAINER"))
    args = parser.parse_args(argv)
    if args.command == "validate-run-config":
        load_config(args.path)
    elif args.command == "validate-locus-config":
        load_locus_config(args.path)
    elif args.command == "validate-case-package":
        validate_case_package(args.path)
    elif args.command == "print-resolved-config":
        print(json.dumps(load_config(args.path), indent=2, sort_keys=True))
    else:
        print(run(args.config,dry_run=args.dry_run,start_stage=args.start_stage,stop_stage=args.stop_stage,resume=args.resume,overwrite=args.overwrite,execution_mode=args.execution_mode))
    return 0
