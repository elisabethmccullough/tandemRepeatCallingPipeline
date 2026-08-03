# tandemRepeatCallingPipeline

A reproducible, **non-clinical** pipeline for collecting tandem-repeat evidence from a local assembly and its original coordinate-sorted mini-BAM.

## Status and planned callers

Stages 03 and 04 implement version-adapted VAMOS read and contig execution. Later callers remain planned. Normalization deliberately does not imply equivalent caller measurements.

## Inputs

A run requires a LoMA-style assembly FASTA, original coordinate-sorted mini-BAM and BAI, reference FASTA and FAI, locus YAML, and filename-safe sample/locus IDs. Inputs are validated as non-empty where appropriate and never modified. Reference assets remain external to Git.

## Install and quick start

Python 3.11 or newer is required.

```bash
python -m pip install -e ".[dev]"
bash scripts/run_pipeline.sh --config config/example.yaml --dry-run
pytest
```

Edit a copy of `config/example.yaml` first. Relative input paths, `locus_config`,
and `run.output_root` are resolved from the directory containing the supplied run
configuration YAML, not from the invocation working directory. Absolute paths are
preserved. This makes a configuration behave identically wherever
`run_pipeline.sh` is invoked. `--start-stage` and `--stop-stage` accept the
numbered names listed below.

## Configuration

```yaml
run:
  sample_id: HG00438
  locus_id: HTT
  output_root: ../outputs
inputs:
  assembly_fasta: path/to/assembly.fasta
  mini_bam: path/to/sample.mini.bam
  mini_bam_index: path/to/sample.mini.bam.bai
  reference_fasta: path/to/hg38.fasta
  reference_fasta_index: path/to/hg38.fasta.fai
locus_config: loci/htt_hg38.yaml
execution:
  threads: 4
  overwrite: false
```

The preliminary HTT file intentionally has null coordinates. Coordinates and tool parameters must be confirmed before real execution.

## Stages and evidence

1. `00_validate_inputs` validates configuration, inputs, identifiers, and reports missing tools as scaffold warnings.
2. `01_prepare_bam` preserves the **raw-read evidence** via links and represents quickcheck/sort/index operations.
3. `02_align_assembly` represents minimap2 alignment and sorted BAM creation from **assembly-derived evidence**.
4. `03_run_vamos_read` uses the original reads; `04_run_vamos_contig` uses the assembled contig alignment.
5. `05_run_straglr` uses original raw reads.
6. `06_prepare_tandem_genotypes` and `07_run_tandem_genotypes` prepare and call from their required alignments.
7. `08_normalize_outputs` maps future native results without manufacturing absent fields.

Raw-read evidence reflects individual alignments and read support. Assembly-derived evidence reflects the assembled sequence and assembly/alignment process. They are complementary, not interchangeable, and provenance remains explicit in `analysis_source`.

## Repository and outputs

Source modules live in `src/tr_calling_pipeline`, entry points in `scripts`, YAML under `config`, tests under `tests`, resources under `resources`, and the future Apptainer recipe under `containers/apptainer`.

Each run creates `outputs/<sample_id>_<locus_id>/` with `00_manifest`, `01_inputs`, `02_prepared_bam`, `03_assembly_alignment`, `04_vamos_read`, `05_vamos_contig`, `06_straglr`, `07_tandem_genotypes`, `08_normalized`, and `logs`. Planned native names are `vamos_read.vcf.gz`, `vamos_contig.vcf.gz`, `straglr.tsv`, `straglr.bed`, and `tandem_genotypes.txt`. Normalized results are `caller_results.tsv`, `caller_results.json`, and `normalization_warnings.tsv`; generated outputs are ignored by Git.

## Reproducibility and limitations

Manifests support resolved paths, SHA-256 checksums, commit identity, tool versions, stage states, outputs, warnings, and completion state. Future work must pin/install caller versions, confirm locus coordinates and caller arguments, implement tool-version capture and lifecycle updates in the orchestrator, implement native parsers, and validate with controlled datasets. Dry-run still requires real non-empty fixture inputs so validation remains meaningful.

This software is for pipeline development and research only. It performs no clinical interpretation and does not classify alleles or provide medical conclusions.

## Versioned pipeline–GUI contracts

JSON Schema Draft 2020-12 contracts live in `schemas/` for run and locus
configuration, lossless normalized caller evidence, and GUI-facing case
manifests. Stable enum values and GUI coordination decisions are documented in
`docs/gui-handoff/task-01-contracts.md`; frozen valid and invalid GUI fixtures
live under `tests/fixtures/case-packages/`.

After installation, validate or inspect contracts with:

```bash
tr-pipeline validate-run-config config/example.yaml
tr-pipeline validate-locus-config config/loci/htt_hg38.yaml
tr-pipeline validate-case-package tests/fixtures/case-packages/minimal
tr-pipeline print-resolved-config config/example.yaml
```

Case-package inventory paths are package-relative POSIX paths. Validation rejects
absolute paths, traversal, duplicate identities, unknown file references, and
symlink escapes. Caller evidence retains raw strings and structured values and
keeps `RAW_READS` separate from `ASSEMBLED_CONTIG`.

## Execution and provenance framework

External processes are represented by immutable argument-vector command specifications and executed without a shell. Tool discovery preserves configured and resolved executable identities, distinguishes required and optional tools, probes raw/normalized versions, and supports stable `NATIVE` and `APPTAINER` modes. Apptainer images, internal commands, and explicit bind mounts are retained in provenance; broad implicit mounts are never added.

Each command streams separate stdout/stderr logs plus a combined convenience log under `logs/<stage>/<command>.*.log`. UTC timestamps, monotonic durations, redacted environment overrides, checksums, failures, timeouts, and output validation are atomically recorded. Existing outputs are rejected unless overwrite is explicit. Required outputs must exist and be non-empty unless declared otherwise. Timeout handling terminates the process group and retains partial diagnostics.

The canonical eleven-stage registry drives selection and stage records. Resume is deliberately conservative: prior success, supported schemas, configuration/input/output checksums, tool/version/mode, and container identity must all match. Dry runs write records and logs but execute nothing. Configuration digests use canonical JSON rather than YAML formatting; stage digests include run/input and relevant tool/container configuration. Configurations without the optional `tools` and `container` sections remain valid.

The scaffold runner implements stage selection and records the observable effect of every public planning option. With `--resume`, a matching successful canonical stage record is preserved and a separate `<stage>.skip.json` record reports `SKIPPED` and `RESUME_ALLOWED`; stale records are archived as `INVALIDATED` before replanning. With `--no-resume`, an existing record is a conflict unless `--overwrite` (or the run configuration's overwrite policy) permits replacement, and the replacement explicitly records `RESUME_DISABLED`. `--execution-mode` overrides configured tool modes in planning provenance. `APPTAINER` requires both a configured launcher and an existing image, whose checksum is recorded. These are scaffold-planning behaviors; the reusable execution module implements process execution, but caller-specific stage commands remain deferred.

```bash
PYTHONPATH=src python -m tr_calling_pipeline.cli validate-run-config config/example.yaml
PYTHONPATH=src python -m tr_calling_pipeline.cli run --config config/example.yaml --dry-run
tr-pipeline run --config config/example.yaml --dry-run
tr-pipeline run --config config/example.yaml --resume
tr-pipeline run --config config/example.yaml --start-stage 02_align_assembly --stop-stage 05_run_straglr
```

Caller-specific stage execution remains incomplete; the current runner plans the scaffold and emits lifecycle provenance without inventing biological parameters.
# Task 3 inputs, BAM preparation, and assembly alignment

The first three stages now validate immutable configured inputs, prepare a stage-local
coordinate-sorted mini-BAM, and align configured patient assembly records. Inputs are
explicit: assembly FASTA records, BAM and its configured index, reference FASTA and
FAI, and a locus configuration. No nearby file is discovered by filename.

`inputs.assembly_records` accepts one or more explicitly selected
`PATIENT_HAPLOTYPE` and/or `PATIENT_CONSENSUS` records. It does not require two
haplotypes and labels such as HAP1/HAP2 do not imply parent or caller allele identity.
FASTA uses the documented IUPAC DNA alphabet `ACGTRYSWKMBDHVN`; letters are retained
exactly, while the biological sequence SHA-256 hashes uppercase, unwrapped letters.
Reference FASTA names and lengths must exactly agree with the configured FAI.
`inputs.reference_scope` may explicitly be `WHOLE_GENOME_REFERENCE`,
`LOCAL_LOCUS_REFERENCE`, or `UNKNOWN_REFERENCE_SCOPE` (the backward-compatible
default).

Stage 01 runs samtools quickcheck, header, flagstat, and idxstats checks. A valid
coordinate-sorted source and index are hard-linked when possible and copied otherwise;
other sort orders are sorted and newly indexed. Sources are never changed. Stage 02
uses separate minimap2, samtools view, sort, and index commands, each with execution
provenance. Its default `asm20` preset is a **development default pending biological
validation**, configurable through `assembly_alignment` with `secondary` and `threads`.

The explicitly configured BAM and index are validated as a pair in a private stage
directory. They are linked or copied to `configured-source.bam` and
`configured-source.bam.bai`, and samtools operates only on that isolated pair; therefore
an unrelated conventionally named index beside the original BAM cannot be selected.
Prepared-BAM `idxstats` likewise operates on `prepared.mini.bam` beside the exact
`prepared.mini.bam.bai`. Existing final prepared outputs are rejected unless overwrite
is explicit. Hard-link fallback is limited to cross-filesystem or unsupported-link
errors; unexpected storage errors remain failures rather than being hidden by a copy.

The authoritative `patient-sequences.fasta` contains selected sequence letters exactly
as stored (apart from deterministic 80-column wrapping) under stable `sequence_id`
headers. Alignment never rewrites or reverse-complements it. Reverse mappings set
`reverse_complement_required: true`; a consumer may transform dynamically for a
reference-forward display. Reverse strand does not imply parental origin. Unresolved
records remain visible with explicit mapping states: `UNIQUE_PRIMARY`,
`MULTIPLE_PRIMARY`, `UNMAPPED`, `AMBIGUOUS`, `TRUNCATED`, `SECONDARY_ONLY`, and
`SUPPLEMENTARY_ONLY`; pre-alignment metadata uses `NOT_ALIGNED`.

```bash
tr-pipeline validate-run-config config/example.yaml
tr-pipeline run --config config/example.yaml --start-stage 00_validate_inputs --stop-stage 02_align_assembly
tr-pipeline run --config config/example.yaml --dry-run --start-stage 00_validate_inputs --stop-stage 02_align_assembly
```

Resume identities include source/output checksums, alignment configuration, and tool
versions. Dry runs create planning/lifecycle records but no BAM or FASTA outputs.
Motif parsing, anchor projection, caller execution, caller allele assignment, clinical
metrics, and clinical interpretation remain intentionally outside this work.

## VAMOS execution (stages 03–04)

VAMOS is not installed automatically. Configure its executable and whether it is required under `tools.vamos`. Catalogs are explicit in the locus configuration (`repeat_catalog`, or separate `read_catalog` and `contig_catalog`) and resolve relative to that file; the pipeline never searches for a substitute. Missing optional tools produce `TOOL_MISSING` metadata; missing required tools fail. Missing catalogs and unsupported versions have explicit statuses.

No real VAMOS release is currently declared supported. The isolated provisional adapter is disabled by default; setting `tools.vamos.allow_provisional_adapter: true` is a development-only opt-in for fake-tool testing with 2.x-shaped version strings. Its syntax still requires verification against the laboratory's installed VAMOS. Every declared native output remains unchanged with checksum, size, caller version, and command provenance. The known JSONL fixture format is normalized losslessly; unrecognized formats remain available and report `UNSUPPORTED_FORMAT`.

Read-derived alleles remain `UNASSIGNED`; allele order is never mapped to patient haplotypes. Assembly records are directly associated with the exact source sequence. VAMOS does not replace the patient FASTA, results are not averaged, and no clinical interpretation is performed. Resume compares input, output, configuration, catalog, and tool-version identities conservatively.

```bash
tr-pipeline run --config config/example.yaml --start-stage 03_run_vamos_read --stop-stage 04_run_vamos_contig
tr-pipeline run --config config/example.yaml --dry-run --start-stage 03_run_vamos_read --stop-stage 04_run_vamos_contig
```

## STRaglr stage (development contract)

Stage `05_run_straglr` consumes the prepared coordinate-sorted mini-BAM and index,
configured reference FASTA and index, and the explicit `caller_resources.straglr.repeat_catalog`
path from the locus YAML. It never searches for or creates a catalog. Configure the
executable under `tools.straglr`; `required: false` records missing tools without
failing the pipeline, while `required: true` fails after preserving metadata.

STRaglr has one canonical configuration split: executable identity, requiredness,
and execution mode are under `tools.straglr`, while the catalog, adapter opt-in,
and caller argv extensions are under
`caller_resources.straglr.{repeat_catalog,allow_provisional_adapter,additional_arguments}`
in the locus configuration. The run-config schema rejects caller-resource settings
under `tools.straglr`, so no accepted value is silently ignored.

No real STRaglr release is currently laboratory verified. The isolated provisional
1.x adapter is disabled by default and requires the locus setting
`caller_resources.straglr.allow_provisional_adapter: true`;
unknown, undetermined, and unsupported versions are never executed with guessed
flags. `additional_arguments` is an argv string list. Its positional syntax and
frozen synthetic TSV layout still require verification against the laboratory
installation.

STRaglr executes into a unique private work directory. Only files from that
execution are atomically published to `06_straglr/native/` and registered with
SHA-256, size, caller version, and command provenance. Without overwrite an existing
completed native set is a conflict; with overwrite it is replaced rather than merged.
Normalization preserves every TSV value as text in `raw_fields`. Evidence uses
`RAW_READS`, remains `UNASSIGNED`, and is never linked to HAP1/HAP2. STRaglr does
not replace patient sequences, values are not reconciled with VAMOS, and no clinical
interpretation is performed. Unsupported native layouts remain registered while
normalization reports `UNSUPPORTED_FORMAT`. Resume includes input, catalog, tool
version, adapter setting, additional arguments, and native output identities.

```bash
tr-pipeline run --config config/example.yaml --start-stage 05_run_straglr --stop-stage 05_run_straglr
tr-pipeline run --config config/example.yaml --dry-run --start-stage 05_run_straglr --stop-stage 05_run_straglr
```

## LAST and tandem-genotypes (development gated)

Stages 06–07 resolve `lastdb`, `lastal`, and `tandem-genotypes` from `tools`.
They may be optional; a missing optional tool yields terminal missing-tool metadata,
whereas required tools fail their stage. Configure the sole repeat resource under
`caller_resources.tandem_genotypes.repeat_definition` in the locus YAML (relative
to that YAML). The pipeline never searches for or invents this file.

There are no laboratory-verified adapters yet. Provisional LAST and caller syntax
must be enabled separately with `allow_provisional_last_adapter` and
`allow_provisional_tandem_genotypes_adapter`; additional argument arrays are part
of resume identity. Unknown versions and formats are preserved and reported, not
guessed. Preparation writes one deterministic FASTA, private LAST database, and
MAF alignment per exact patient sequence. Native caller files remain immutable;
normalization retains raw strings and direct sequence association. Unverified
coordinates remain null with `UNKNOWN_COORDINATE_SPACE`.

```bash
tr-pipeline run --config config/example.yaml --start-stage 06_prepare_tandem_genotypes --stop-stage 07_run_tandem_genotypes
tr-pipeline run --config config/example.yaml --dry-run --start-stage 06_prepare_tandem_genotypes --stop-stage 07_run_tandem_genotypes
```

LAST alignment never modifies the authoritative sequence. No parental origin or
clinical interpretation is inferred. Real integration remains pending tests on
verified laboratory installations.
