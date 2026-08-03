# Final project summary

## What it does
The non-clinical pipeline accepts a patient assembly FASTA (commonly produced upstream by LoMA), mini-BAM/index, reference/index, locus configuration, and caller resources. It validates immutable inputs, prepares reads, aligns patient sequences, invokes explicitly configured callers, preserves native files, normalizes distinct evidence, constructs a portable GUI package, and validates that package independently.

## Outputs
Outputs include patient sequences, separate caller evidence, native files, provenance, warnings, checksums, and a portable GUI package. The GUI consumes `case-manifest.json` and referenced relative files, never arbitrary run-directory state.

## Current verification and remaining work
Internal contracts and adapters have unit/fake-tool synthetic coverage. Fake executables are not real callers. No real caller installation or laboratory workflow is verified. Remaining work is controlled real-tool smoke testing, laboratory review, final HTT resource confirmation, approved real-sample evaluation, GUI integration, and clinical governance if ever pursued.

## Scientific limitations
There is no interpretation, thresholding, consensus, caller preference, parent-of-origin inference, or raw-read haplotype assignment. The package is technical evidence, not a clinical report. Current status: **READY_WITH_LIMITATIONS** for continued development and synthetic handoff, not clinical use.
