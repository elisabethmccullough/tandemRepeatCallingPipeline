# GUI handoff contract index

| Contract | Documentation | Fixture | GUI entry file |
|---|---|---|---|
| Patient sequences | `task-03-sequences-bam-and-orientation.md` | `tests/fixtures/gui-handoff/task-03` | patient metadata referenced by manifest |
| VAMOS | `task-04-vamos-read-and-assembly.md` | `task-04` | normalized evidence |
| STRaglr | `task-05-straglr.md` | `task-05` | normalized evidence |
| tandem-genotypes | `task-06-tandem-genotypes.md` | `task-06` | normalized evidence |
| Unified normalization | `task-07-unified-normalization.md` | `task-07` | `normalized-evidence.json` |
| Portable package | `task-08-case-package.md` | `tests/fixtures/case-packages` | `case-manifest.json` |

Version 1.x schema fields and enum meanings are stable only as documented. Adapter command syntax remains provisional. The GUI may display evidence, sequence association explicitly present in records, provenance, fake/real executable kind, verification state, and warnings. It must not invent consensus, clinical classifications, caller preference, parental origin, read-to-sequence assignment, or rewrite native data. Validate the package before opening it.
