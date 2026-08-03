# Task 07: unified normalization GUI handoff

## Purpose

The unified normalization package gives the GUI one consistent place to load
all caller evidence. It standardizes field names while preserving caller
identity and version, raw fields, native-file provenance, exact patient
sequence associations, warnings, and missing or unsupported states.

Normalization **does not** make a consensus, choose a preferred caller,
average values, assign read calls to patient haplotypes, classify expansions,
or make clinical interpretations. Shared fields are common display slots, not
proof that callers measured equivalent biological quantities.

## GUI grouping guidance

The GUI may group records as **Unassigned read evidence**, **Patient sequence
1**, and **Patient sequence 2**. Caller-specific records remain separate inside
each group. Sequence labels are stable identifiers and do not imply parental
origin.

| Caller mode | Analysis source | Assignment |
|---|---|---|
| VAMOS read | `RAW_READS` | `UNASSIGNED` |
| STRaglr | `RAW_READS` | `UNASSIGNED` |
| VAMOS contig | `ASSEMBLED_CONTIG` | `DIRECT_SEQUENCE_ASSOCIATION` |
| tandem-genotypes | `ASSEMBLED_CONTIG` | `DIRECT_SEQUENCE_ASSOCIATION` |

Records sort by analysis source, sequence ID (null first), caller, native record
ID, native allele ID, and record ID. Sources sort by caller, analysis source,
sequence ID, and file ID. Caller summaries use VAMOS-read, STRaglr,
VAMOS-contig, tandem-genotypes order.

The GUI should display, or make accessible, unsupported-format and
missing/failed-caller states, provisional-adapter and unverified-coordinate
warnings, and source/version provenance. Availability summaries are inventory
only and never claim agreement or disagreement.

