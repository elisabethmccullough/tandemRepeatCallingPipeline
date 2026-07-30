#!/usr/bin/env bash
set -euo pipefail
# Preserve the raw-read mini-BAM by linking it; future releases will verify/sort only when needed.
usage(){ echo "Usage: $0 --bam FILE --bai FILE --output-dir DIR [--samtools CMD] [--dry-run]"; }
bam=""; bai=""; out=""; samtools=samtools; dry=false
while (($#)); do case "$1" in --bam) bam=${2:-}; shift 2;; --bai) bai=${2:-}; shift 2;; --output-dir) out=${2:-}; shift 2;; --samtools) samtools=${2:-}; shift 2;; --dry-run) dry=true; shift;; --help|-h) usage; exit 0;; *) echo "Error: unknown argument: $1" >&2; exit 2;; esac; done
for value in bam bai out; do [[ -n ${!value} ]] || { echo "Error: --${value/output-dir} is required" >&2; exit 2; }; done
[[ -s "$bam" && -f "$bai" ]] || { echo "Error: BAM must be non-empty and index must exist" >&2; exit 1; }
mkdir -p "$out"; log="$out/prepare_bam.log"
commands=("$samtools quickcheck '$bam'" "ln -s '$(realpath "$bam")' '$out/original.mini.bam'" "ln -s '$(realpath "$bai")' '$out/original.mini.bam.bai'")
printf 'Intended commands (sorting, if required: %s sort; indexing: %s index):\n' "$samtools" "$samtools" | tee "$log"
for cmd in "${commands[@]}"; do echo "+ $cmd" | tee -a "$log"; $dry || eval "$cmd"; done
