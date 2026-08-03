"""Command-line interface for shared contract validation."""

from __future__ import annotations

import argparse
import json
from typing import Sequence

from .case_package_validation import validate_case_package
from .config import load_config, load_locus_config
from .runner import run


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="tr-pipeline")
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("validate-run-config", "validate-locus-config", "print-resolved-config"):
        command = commands.add_parser(name)
        command.add_argument("path")
    command = commands.add_parser("validate-case-package")
    command.add_argument("--package", required=True)
    command = commands.add_parser("run")
    command.add_argument("--config", required=True)
    command.add_argument("--dry-run", action="store_true")
    command.add_argument("--start-stage")
    command.add_argument("--stop-stage")
    command.add_argument("--resume", action=argparse.BooleanOptionalAction, default=True)
    command.add_argument("--overwrite", action="store_true")
    command.add_argument("--execution-mode", choices=("NATIVE", "APPTAINER"))
    for name in ("validate-inputs", "prepare-bam", "align-assembly", "build-case-package"):
        focused = commands.add_parser(name)
        focused.add_argument("--config", required=True)
        focused.add_argument("--dry-run", action="store_true")
        focused.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    if args.command == "validate-run-config":
        load_config(args.path)
    elif args.command == "validate-locus-config":
        load_locus_config(args.path)
    elif args.command == "validate-case-package":
        report = validate_case_package(args.package, write_report=True)
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0 if report["valid"] else 2
    elif args.command == "print-resolved-config":
        print(json.dumps(load_config(args.path), indent=2, sort_keys=True))
    elif args.command == "run":
        print(run(args.config,dry_run=args.dry_run,start_stage=args.start_stage,stop_stage=args.stop_stage,resume=args.resume,overwrite=args.overwrite,execution_mode=args.execution_mode))
    else:
        stage = {"validate-inputs":"00_validate_inputs", "prepare-bam":"01_prepare_bam", "align-assembly":"02_align_assembly", "build-case-package":"09_build_case_package"}[args.command]
        print(run(args.config, dry_run=args.dry_run, start_stage=stage, stop_stage=stage, resume=True, overwrite=args.overwrite))
    return 0
