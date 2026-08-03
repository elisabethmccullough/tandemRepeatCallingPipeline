# Task 05: STRaglr GUI handoff (contract 1.0)

## Stable caller identity

`caller: STRAGLR` is stable.

## Stable evidence distinction

STRaglr evidence is read-derived and uses `analysis_source: RAW_READS`. Calls remain
`UNASSIGNED`; `associated_sequence_id` is null. STRaglr alleles are not patient
haplotypes, native order must not be mapped to HAP1/HAP2, and STRaglr does not
replace the patient FASTA.

## Stable fields

`caller`, `caller_version`, `analysis_source`, `native_record_identifier`,
`native_allele_identifier`, `associated_sequence_id`, `assignment_state`,
`source_file_id`, `source_file_sha256`, `raw_fields`, and
`normalization_warnings` are stable serialized fields.

## GUI display guidance

The GUI may group STRaglr as unassigned read evidence, display native allele IDs,
caller version and exact native values, show warnings, and link to native files.
It must not present allele 1/2 as HAP1/HAP2, infer parental origin, average STRaglr
with VAMOS, hide native values, or present clinical interpretation.

## VAMOS versus STRaglr

Both may provide read-derived evidence, but can report different native data types.
Their values and native evidence roles remain separate. No equivalence,
prioritization, or reconciliation is implied; comparison is future work.

## Provisional decisions

No STRaglr version is laboratory verified. The isolated 1.x development adapter,
its positional command spelling, frozen TSV column layout, catalog format, and
zero-based half-open coordinate convention are provisional and require explicit
opt-in. Catalog contents are opaque to the pipeline. Future allele assignment,
display columns, and cross-caller reconciliation are deferred.

## Matching GUI work

The GUI can begin an unassigned STRaglr read-evidence group, native-field
inspection, caller-version display, unsupported-format warnings, and separate
VAMOS/STRaglr evidence roles using the frozen fixtures.
