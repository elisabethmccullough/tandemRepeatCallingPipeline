# Task 02: execution and provenance handoff

## Stable values

Tool IDs are `SAMTOOLS`, `MINIMAP2`, `VAMOS`, `STRAGLR`, `LASTDB`, `LASTAL`, `TANDEM_GENOTYPES`, `PYTHON`, `PIPELINE`, and `APPTAINER`. Tool statuses are `AVAILABLE`, `MISSING_OPTIONAL`, `MISSING_REQUIRED`, `VERSION_UNDETERMINED`, `UNSUPPORTED_VERSION`, and `NOT_CHECKED`. Execution modes are `NATIVE` and `APPTAINER`. Execution statuses are `PLANNED`, `RUNNING`, `SUCCEEDED`, `FAILED`, `TIMED_OUT`, `CANCELLED`, `DRY_RUN`, and `OUTPUT_VALIDATION_FAILED`. Stage statuses are `NOT_STARTED`, `PLANNED`, `RUNNING`, `SUCCEEDED`, `FAILED`, `SKIPPED`, `DRY_RUN`, and `INVALIDATED`.

Resume reasons are `RESUME_ALLOWED`, `RESUME_DISABLED`, `NO_PRIOR_RECORD`, `PRIOR_STAGE_FAILED`, `INPUT_CHANGED`, `OUTPUT_MISSING`, `OUTPUT_CHANGED`, `CONFIGURATION_CHANGED`, `TOOL_CHANGED`, `TOOL_VERSION_CHANGED`, `EXECUTION_MODE_CHANGED`, `CONTAINER_CHANGED`, and `UNSUPPORTED_RECORD_VERSION`. Tool-status, execution-record, and stage-record schemas are version `1.0`. Enum values serialize as the strings above, never language-specific enum representations. The optional run-schema `tools` and `container` additions are pipeline-only and backward compatible at run schema `1.0`.

## Pipeline artifacts available to the GUI

Frozen fixtures expose pipeline version, optional Git commit and dirty state, tool/caller display names and versions, complete input/output SHA-256 identities, execution mode and container identity, stage status, and structured failure status. Raw version output remains diagnostic; normalized caller versions are the display-oriented value. These records supplement rather than overload the case manifest.

## Visible PDF provenance recommendation

Visible provenance may include pipeline name/version, Git commit when available, caller names/versions, execution mode when relevant, abbreviated source checksums, and case-package creation time. Full commands and local paths should not be required or shown.

## Embedded diagnostic provenance recommendation

Machine-readable provenance may embed full argument vectors, working directories, configured/resolved tool paths, container image hashes, bind mounts, complete checksums, stage/execution records, bounded failure diagnostics, and configuration digests.

## Privacy and portability warning

Absolute paths are workstation-specific and must not be displayed in clinical reports. The GUI should sanitize or omit them in visible PDF content and prefer portable package-relative file identities. Detailed records may retain paths for diagnostics and therefore require appropriate handling.

## Provisional decisions

Caller-specific argument layouts, supported-version ranges, clinical relevance of execution mode, stage-to-package inventory placement, and final caller-stage output roles remain provisional. The GUI must not treat raw local paths, full command rendering, bind layout, or scaffold warnings as stable presentation fields. No evidence classification or caller precedence is introduced.

The scaffold resume policy preserves a canonical successful stage record and emits a sibling `<stage>.skip.json` record for each allowed skip. Incompatible canonical records are archived with `INVALIDATED` status before a new plan is written. `RESUME_DISABLED` explicitly identifies planning performed with `--no-resume`; it is stable as a resume reason in record schema version `1.0`.

## Matching GUI work

The GUI can begin schema/fixture validation, pipeline and caller version display, checksum identity display, failure diagnostics, provenance embedding, and visible-PDF provenance design using `tests/fixtures/provenance/`. Real caller rendering should wait for caller-specific stages and normalized evidence integration.
