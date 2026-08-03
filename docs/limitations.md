# Limitations and scientific boundary

This pipeline performs no clinical interpretation, caller consensus, averaging, clinical thresholds, parent-of-origin inference, or automatic raw-read-to-haplotype assignment. Caller evidence remains separate. Adapters and native formats are provisional; no real versions or laboratory workflows have been verified, and fixtures are synthetic only. HTT resources may remain unresolved or provisional, and laboratory caller formats may differ.

The case package is a technical evidence package, not a clinical report. LoMA is upstream and file-based; it is not directly run unless a future explicit integration is added:

`Long reads → LoMA or another assembler → patient sequence FASTA → tandem-repeat caller pipeline → portable GUI case package`

Real-sample work requires approved processes, resource confirmation, tool verification, and appropriate governance.
