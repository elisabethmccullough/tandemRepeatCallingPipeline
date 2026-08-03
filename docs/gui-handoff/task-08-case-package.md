# Task 8 GUI handoff: portable case package

The case package is the complete portable folder opened by the GUI. Before
opening it as trusted, the GUI inspects `case-manifest.json` and
`package-validation.json`; an unsuccessful report must produce a visible
integrity warning identifying changed or missing files.

## Authoritative patient data

`patient/patient-sequences.fasta` and
`patient/patient-sequences.metadata.json` are authoritative. Caller output
never replaces these upstream-assembler (for example, LoMA) sequences. The GUI
loads unified evidence from `evidence/normalized-evidence.json`,
`evidence/evidence-summary.json`, and `evidence/source-registry.json`. It may
link to `native/`, but normalized and native values remain separate.

## Portability and integrity

All GUI-required paths are POSIX-style paths relative to the package root.
Original source paths are provenance only. The checksum inventory covers every
payload plus the manifest; it explicitly excludes itself and the regenerable
validation report. Extra files are validation errors. Regular files are copied
rather than linked, and traversal, symlinks, special files, duplicate paths,
checksum changes, missing files, and secret-like configuration keys are
rejected.

## Stable limitations

The package has no consensus call, preferred caller, parental-origin inference,
expansion classification, clinical threshold, or clinical interpretation. It
can retain unavailable, failed, unsupported, and provisional-adapter evidence
and warnings without treating them as biological conclusions.

