# Operational runbook

## Before running
Confirm the assembly FASTA, mini-BAM and exact index, reference and index, locus configuration, explicit caller catalogs/definitions, tool versions, and output permissions. Record checksums; do not edit authoritative inputs.

## Dry run and full run
```bash
tr-pipeline run --config "/path with spaces/run.yaml" --dry-run
tr-pipeline run --config "/path with spaces/run.yaml"
```
Dry-run plans without biological execution. A full run requires explicitly installed tools; nothing is downloaded.

## Resume, overwrite, and optional callers
Resume is conservative and accepts only matching configuration, inputs, outputs, and tool identities. Changed upstream identities invalidate affected downstream work. `--no-resume --overwrite` explicitly replaces conflicting stage outputs; archive a valid package first. Missing optional callers terminate in an explicit unavailable state and must not generate evidence.

## Troubleshooting
- **Tool missing/version undetermined/unsupported adapter:** correct the explicit executable or version; never guess flags. Provisional adapters require opt-in.
- **Missing catalog/definition:** configure the exact local resource; no discovery occurs.
- **Checksum/native-output conflict:** preserve files, investigate identity changes, then use explicit overwrite only when intended.
- **Package validation failure:** inspect issue codes; do not hand an invalid directory to the GUI.
- **Spaces/Windows/WSL:** quote shell paths. Source paths may use platform syntax, but package JSON paths are relative POSIX paths. WSL is not a separately supported CI platform.
- **Apptainer unavailable:** use verified native execution or install/configure it outside the pipeline; the pipeline never installs it.

## GUI handoff
Provide the completed case-package root containing `case-manifest.json` only after `validate-case-package` returns valid. The original run directory is not the GUI contract.
