# Task 03: sequences, BAM preparation, and orientation handoff

## Stable sequence fields

Version 1.0 freezes `sequence_id`, `sequence_role`, `display_label`,
`source_fasta_record_id`, `sequence_sha256`, `sequence_length`,
`original_orientation`, `reference_alignment_strand`, `display_orientation`,
`reverse_complement_required`, `source_coordinates`, and `mapping_status`.
`sequence_id` is the package FASTA record ID; the source ID remains explicit provenance.

## Stable enum values

* `original_orientation`: `AS_STORED`.
* `reference_alignment_strand`: `FORWARD`, `REVERSE`, `UNKNOWN`, `NOT_ALIGNED`.
* `display_orientation`: `REFERENCE_FORWARD`, `UNRESOLVED`.
* `mapping_status`: `UNIQUE_PRIMARY`, `MULTIPLE_PRIMARY`, `UNMAPPED`,
  `AMBIGUOUS`, `TRUNCATED`, `SECONDARY_ONLY`, `SUPPLEMENTARY_ONLY`, `NOT_ALIGNED`.

A unique primary has exactly one non-secondary, non-supplementary mapped record.
Multiple primary, secondary-only, and supplementary-only have their literal SAM-derived
meanings. Unmapped has no mapped record. Ambiguous has mapped evidence that does not
fit a single category. `TRUNCATED` is reserved for later validated completeness checks;
`NOT_ALIGNED` is the explicit pre-alignment state. Coordinates are always
`zero_based_half_open` (SAM POS minus one through start plus reference-consuming CIGAR
length).

## Authoritative FASTA policy

`03_assembly_alignment/patient-sequences.fasta` is authoritative. It retains exact
selected as-stored letters and deterministic wrapping; it is never reverse-complemented.
The sequence hash is SHA-256 of uppercase, unwrapped IUPAC letters, so case and wrapping
do not change biological identity. The original source-file SHA-256 independently
preserves byte identity. The GUI should dynamically reverse-complement only for a
forward-locus view when `reverse_complement_required` is true. A future derived display
FASTA, if adopted, must have a separate identity and provenance.

## GUI responsibilities

The GUI may load the package FASTA, match stable sequence IDs, retain source IDs,
reverse-complement for display, show mapping warnings, and keep unmapped sequences
visible. It must not infer parental origin, interpret HAP1/HAP2 as maternal/paternal,
replace exact sequence with caller output, silently recalculate mapping status, or hide
unmapped/ambiguous evidence without warning.

## Provisional decisions

* Minimap2 `asm20` is only a development default pending biological validation.
* Dynamic display transformation currently belongs to the GUI; ownership may become a
  future shared derived-artifact contract.
* Final case-manifest references are deferred to package construction.
* Anchor-derived orientation checks and cross-build transformations are future work.

## Matching GUI work

The frozen fixtures under `tests/fixtures/gui-handoff/task-03` can now drive loading of
patient sequences, orientation metadata, mappings, exact checksums, and warnings. They
include synthetic forward, reverse, consensus/not-aligned, and unmapped examples and no
patient data.
