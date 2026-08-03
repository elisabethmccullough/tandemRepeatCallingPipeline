# Task 6 tandem-genotypes handoff

## Stable identity and evidence distinction

`caller` is `TANDEM_GENOTYPES`; `analysis_source` is `ASSEMBLED_CONTIG`; and
`assignment_state` is `DIRECT_SEQUENCE_ASSOCIATION`. Every record is associated
with the exact configured patient sequence used for LAST preparation and caller
execution. The caller output does not replace the authoritative patient FASTA,
and LAST alignment is an intermediate evidence artifact. Sequence IDs do not
imply parental origin.

Stable evidence fields are `caller`, `caller_version`, `analysis_source`,
`native_record_identifier`, `native_allele_identifier`, `associated_sequence_id`,
`assignment_state`, `source_file_id`, `source_file_sha256`, `raw_fields`, and
`normalization_warnings`. Alignment metadata exposes `sequence_id`,
`source_fasta_record_id`, `sequence_sha256`, `alignment_file_id`,
`alignment_file_sha256`, `alignment_format`, `coordinate_space`, `lastdb_version`,
`lastal_version`, and `repeat_definition_identity`.

## GUI guidance

The GUI may display records below their associated sequence, native motif and
allele fields, alignment provenance, versions, warnings, and native-file links.
It must not treat caller output as exact sequence, infer parents, merge callers,
hide native values, offer clinical interpretation, or label alignment coordinates
as reference coordinates without explicit metadata.

VAMOS read mode and STRaglr are unassigned raw-read evidence. VAMOS contig mode
and tandem-genotypes are directly associated assembly evidence; tandem-genotypes
follows LAST preparation. This structure implies no equivalence or priority.

## Provisional decisions

There are currently no laboratory-verified LAST or tandem-genotypes versions.
Numeric LAST releases and tandem-genotypes 0.x use development-only adapters,
disabled by default. The frozen formats are MAF and development TSV v1.
Repeat-definition syntax and coordinate semantics require laboratory verification;
unknown coordinates remain null and use `UNKNOWN_COORDINATE_SPACE`. Cross-caller
comparison is deferred.

The GUI can begin sequence-bound tracks, LAST provenance, native-field and version
inspection, unsupported-format warnings, and separate assembly-evidence roles.
