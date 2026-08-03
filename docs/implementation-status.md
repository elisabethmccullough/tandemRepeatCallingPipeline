# Implementation and verification status

Verification levels are `UNIT_TESTED`, `SYNTHETIC_INTEGRATION_TESTED`, `REAL_TOOL_SMOKE_TESTED`, `LABORATORY_VERIFIED`, and `UNVERIFIED`. Fake-tool tests never grant a real-tool level.

| Component | Implementation | Verification | Versions / fake coverage / real coverage | Limitation and next action |
|---|---|---|---|---|
| Input validation | Complete | SYNTHETIC_INTEGRATION_TESTED | internal 0.1; fixtures; no external tool | Confirm approved laboratory inputs. |
| Mini-BAM preparation | Complete | SYNTHETIC_INTEGRATION_TESTED | fake samtools; none real | Smoke-test supported samtools. |
| Assembly alignment | Complete | SYNTHETIC_INTEGRATION_TESTED | fake minimap2/samtools; none real | Verify parameters and real layouts. |
| Patient-sequence packaging | Complete | SYNTHETIC_INTEGRATION_TESTED | two-sequence fixtures; n/a | GUI acceptance testing. |
| VAMOS read adapter | Complete, development-gated | SYNTHETIC_INTEGRATION_TESTED | fake 2.x; none real | Syntax/layout require laboratory verification. |
| VAMOS contig adapter | Complete, development-gated | SYNTHETIC_INTEGRATION_TESTED | fake 2.x; none real | Syntax/layout require laboratory verification. |
| STRaglr adapter | Complete, development-gated | SYNTHETIC_INTEGRATION_TESTED | fake 1.x; none real | Syntax/layout require laboratory verification. |
| LASTDB adapter | Complete, development-gated | SYNTHETIC_INTEGRATION_TESTED | fake tool; none real | Verify real version and MAF output. |
| LASTAL adapter | Complete, development-gated | SYNTHETIC_INTEGRATION_TESTED | fake tool; none real | Verify real version and MAF output. |
| tandem-genotypes adapter | Complete, development-gated | SYNTHETIC_INTEGRATION_TESTED | fake tool; none real | Verify command and native layout. |
| VAMOS / STRaglr / tandem-genotypes parsing | Complete for frozen formats | SYNTHETIC_INTEGRATION_TESTED | fake native fixtures; none real | Compare real native output without modifying it. |
| Unified normalization | Complete | SYNTHETIC_INTEGRATION_TESTED | all fake caller groups; n/a | Measurements remain distinct; no consensus. |
| Case-package construction / independent validation | Complete | SYNTHETIC_INTEGRATION_TESTED | moved-package tests; n/a | Technical evidence package only. |
| Native execution | Complete framework | UNIT_TESTED | synthetic commands; none bioinformatics | Real-tool smoke tests required. |
| Apptainer execution | Planning complete | UNIT_TESTED | argv/provenance tests; none real | Environment verification required. |
| Resume / dry-run | Complete, conservative | UNIT_TESTED | pytest coverage; n/a | Extend full-run invalidation coverage. |
| GUI handoff contracts | Complete for v1 contracts | SYNTHETIC_INTEGRATION_TESTED | frozen fixtures; no GUI test | GUI must not infer interpretation/consensus. |
| CI | Implemented; hosted result required | UNIT_TESTED until Actions succeeds | Ubuntu/Windows, Python 3.11/3.12 configured | Workflow configuration alone is not platform verification; inspect PR checks. |
| Real-tool verification | Not completed | UNVERIFIED | verified versions: none | Complete controlled templates and review. |
