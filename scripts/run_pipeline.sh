#!/usr/bin/env bash
set -euo pipefail
# Sequential scaffold orchestrator. External caller stages only describe intended commands.
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repository_root=$(cd -- "$script_dir/.." && pwd)
export PYTHONPATH="$repository_root/src${PYTHONPATH:+:$PYTHONPATH}"
usage(){ echo "Usage: $0 --config FILE [--dry-run] [--start-stage STAGE] [--stop-stage STAGE]"; }
config=""; dry=false; start="00_validate_inputs"; stop="08_normalize_outputs"
while (($#)); do case "$1" in --config) config=${2:-}; shift 2;; --dry-run) dry=true; shift;; --start-stage) start=${2:-}; shift 2;; --stop-stage) stop=${2:-}; shift 2;; --help|-h) usage; exit 0;; *) echo "Error: unknown argument: $1" >&2; usage >&2; exit 2;; esac; done
[[ -n "$config" ]] || { echo "Error: --config is required" >&2; exit 2; }; [[ -f "$config" ]] || { echo "Error: configuration not found: $config" >&2; exit 1; }
config=$(realpath "$config")
mapfile -t values < <(python - "$config" <<'PY'
import sys
from tr_calling_pipeline.config import load_config
c=load_config(sys.argv[1])
for v in (c['run']['sample_id'],c['run']['locus_id'],c['run']['output_root'],c['inputs']['assembly_fasta'],c['inputs']['mini_bam'],c['inputs']['mini_bam_index'],c['inputs']['reference_fasta'],c['locus_config'],str(c['execution']['threads'])): print(v)
PY
)
sample=${values[0]}; locus_id=${values[1]}; root=${values[2]}; assembly=${values[3]}; bam=${values[4]}; bai=${values[5]}; reference=${values[6]}; locus_config=${values[7]}; threads=${values[8]}
run="$root/${sample}_${locus_id}"; mkdir -p "$run"/{00_manifest,01_inputs,02_prepared_bam,03_assembly_alignment,04_vamos_read,05_vamos_contig,06_straglr,07_tandem_genotypes,08_normalized,logs}
# Establish provenance files before any stage runs. Stage lifecycle enrichment is future work.
python - "$config" "$run/00_manifest" <<'PY'
import csv, shutil, sys
from pathlib import Path
from tr_calling_pipeline.config import load_config
from tr_calling_pipeline.manifest import create_manifest, write_manifest
cfg = load_config(sys.argv[1])
directory = Path(sys.argv[2])
manifest = create_manifest(cfg, sys.argv[1])
write_manifest(manifest, directory / "run_manifest.json")
with (directory / "input_checksums.tsv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t"); writer.writerow(["input", "sha256"])
    writer.writerows(sorted(manifest["input_checksums"].items()))
with (directory / "tool_versions.tsv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t"); writer.writerow(["tool", "executable", "status"])
    for name, executable in sorted(cfg.get("tools", {}).items()):
        writer.writerow([name, executable, "available" if shutil.which(executable) else "not_found"])
PY
stages=(00_validate_inputs 01_prepare_bam 02_align_assembly 03_run_vamos_read 04_run_vamos_contig 05_run_straglr 06_prepare_tandem_genotypes 07_run_tandem_genotypes 08_normalize_outputs)
index(){ local needle=$1; for i in "${!stages[@]}"; do [[ ${stages[$i]} == "$needle" ]] && { echo "$i"; return; }; done; echo "Error: unknown stage: $needle" >&2; return 1; }
first=$(index "$start"); last=$(index "$stop"); (( first <= last )) || { echo "Error: start stage follows stop stage" >&2; exit 2; }
dryarg=(); $dry && dryarg=(--dry-run)
commands=(
"'$script_dir/00_validate_inputs.sh' --config '$config' ${dryarg[*]}"
"'$script_dir/01_prepare_bam.sh' --bam '$bam' --bai '$bai' --output-dir '$run/02_prepared_bam' ${dryarg[*]}"
"'$script_dir/02_align_assembly.sh' --assembly '$assembly' --reference '$reference' --output-dir '$run/03_assembly_alignment' --threads '$threads' ${dryarg[*]}"
"'$script_dir/03_run_vamos_read.sh' --input '$bam' --locus-config '$locus_config' --output-dir '$run/04_vamos_read' ${dryarg[*]}"
"'$script_dir/04_run_vamos_contig.sh' --input '$assembly' --locus-config '$locus_config' --output-dir '$run/05_vamos_contig' ${dryarg[*]}"
"'$script_dir/05_run_straglr.sh' --input '$bam' --locus-config '$locus_config' --output-dir '$run/06_straglr' ${dryarg[*]}"
"'$script_dir/06_prepare_tandem_genotypes.sh' --input '$bam' --locus-config '$locus_config' --output-dir '$run/07_tandem_genotypes' ${dryarg[*]}"
"'$script_dir/07_run_tandem_genotypes.sh' --input '$run/07_tandem_genotypes' --locus-config '$locus_config' --output-dir '$run/07_tandem_genotypes' ${dryarg[*]}"
"python '$script_dir/08_normalize_outputs.py' --output-dir '$run/08_normalized' ${dryarg[*]}"
)
for ((i=first;i<=last;i++)); do log="$run/logs/${stages[$i]}.log"; echo "==> ${stages[$i]}"; eval "${commands[$i]}" 2>&1 | tee "$log"; done
echo "Pipeline scaffold completed: $run"
