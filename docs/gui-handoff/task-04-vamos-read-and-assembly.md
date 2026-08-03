# Task 04: VAMOS read and assembly handoff

## Stable caller identity and evidence distinctions

`caller` is `VAMOS`. Read mode is `RAW_READS`; assembly mode is
`ASSEMBLED_CONTIG`. Read alleles remain `UNASSIGNED`. Assembly evidence is
`DIRECT_SEQUENCE_ASSOCIATION` to the exact stable sequence ID on which it ran.
VAMOS output does not replace the patient FASTA. The GUI must show read and
assembly evidence as separate roles.

## Stable fields

`record_id`, `caller`, `caller_version`, `analysis_source`,
`native_record_identifier`, `native_allele_identifier`,
`associated_sequence_id`, `assignment_state`, `source_file_id`,
`source_file_sha256`, `raw_fields`, and `normalization_warnings` are stable 1.0
fields. Native fields and numeric-looking strings are retained losslessly.

## GUI display guidance

The GUI may bind assembly evidence beneath its associated patient sequence,
group read evidence as unassigned, show exact values and provenance, and warn
for unsupported normalization. It must not treat VAMOS as the exact patient
sequence, force allele 1/2 to HAP1/HAP2, hide native values, average evidence,
or present a clinical interpretation.

## Provisional decisions

Only version strings matching VAMOS 2.x select the isolated provisional
adapter. Its command flags, final native formats, and catalog format require
confirmation against the laboratory installation. Read-derived assignment,
caller-specific display columns, and cross-caller comparison remain deferred.

## Matching GUI work

Assembly-track binding, unassigned read grouping, caller-version display,
native-field inspection, and unsupported-format warnings can now use the frozen
fixtures in `tests/fixtures/gui-handoff/task-04`.
