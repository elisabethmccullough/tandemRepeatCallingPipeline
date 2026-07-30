#!/usr/bin/env python3
"""Create the common empty result files until caller-specific parsers exist."""
import argparse
from tr_calling_pipeline.normalize import write_normalized_outputs

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if args.dry_run:
        print(f"Would create normalized result files in {args.output_dir}")
    else:
        write_normalized_outputs(args.output_dir)
        print("Placeholder: wrote empty normalized outputs; caller parsers are not enabled yet.")
if __name__ == "__main__": main()
