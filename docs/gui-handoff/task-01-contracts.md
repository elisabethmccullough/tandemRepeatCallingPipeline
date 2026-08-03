# Task 01: pipeline–GUI contracts

## Frozen GUI inputs

The GUI may begin package validation against JSON Schema Draft 2020-12. Consume
`schemas/case-manifest.schema.json` and
`schemas/normalized-caller-evidence.schema.json`; the pipeline also owns the run
and locus schemas beside them. Contract version `1.0` is serialized in every
schema and document.

Frozen package fixtures are under `tests/fixtures/case-packages/`: `minimal`
contains HAP1, HAP2, and consensus entries; `multi-locus` demonstrates two locus
assessments; and `invalid-traversal` must be rejected. Tests additionally create
absolute-path, unknown-reference, duplicate-ID, and symlink-escape cases.

## Stable vocabularies

* Sequence roles: `PATIENT_HAPLOTYPE`, `PATIENT_CONSENSUS`,
  `REFERENCE_SEQUENCE`, `COMPARISON_SEQUENCE`. Pipeline assembly inputs produce
  only the first two.
* Analysis sources: `RAW_READS`, `ASSEMBLED_CONTIG`.
* Assignment states: `DIRECT_SEQUENCE_ASSOCIATION`, `UNASSIGNED`,
  `ASSIGNMENT_NOT_APPLICABLE`. Raw-read inference is not supported.
* Evidence states: `AVAILABLE`, `NOT_APPLICABLE`, `INPUT_MISSING`,
  `COMPUTATION_FAILED`, `AMBIGUOUS`, `NOT_COMPUTED`, `UNSUPPORTED_FORMAT`.
* Coordinate conventions: `zero_based_half_open`, `one_based_closed`,
  `inter_base`.

## Provisional fields and matching decisions

The contents of each locus assessment permit additional properties so GUI-only
locus fields can evolve without being invented by the pipeline. The
`gui_contract` section of a locus configuration is likewise an explicit
extension point. Comparison-panel content, final caller-specific `raw_fields`,
orientation display policy, and validated biological coordinates/anchors remain
provisional. The GUI must preserve decimal-looking shared values as their exact
serialized strings when strings are supplied, retain structured motif chains,
and must not infer haplotype, parental origin, interpretation, or caller
priority. Package paths are POSIX package-relative paths; both implementations
must reject absolute paths, traversal, duplicate identities, unknown file
references, and filesystem symlink escapes.
