"""Lossless, caller-neutral packaging of preliminary caller evidence.

Shared fields are display slots only.  This module deliberately performs no
consensus calling, numerical comparison, haplotype inference, or clinical
interpretation; caller adapters remain responsible for parsing native data.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json
import os
from pathlib import Path
import shutil
import uuid
from typing import Any, Iterable, Mapping, Sequence

from .provenance import atomic_write_json, canonical_digest, sha256_file, utc_now
from .schema_validation import SchemaViolation, validate

SCHEMA_VERSION = "1.0"
CALLER_ORDER = (("VAMOS", "RAW_READS"), ("STRAGLR", "RAW_READS"),
                ("VAMOS", "ASSEMBLED_CONTIG"), ("TANDEM_GENOTYPES", "ASSEMBLED_CONTIG"))
EVIDENCE_STATES = {"AVAILABLE", "NOT_APPLICABLE", "INPUT_MISSING", "COMPUTATION_FAILED",
                   "AMBIGUOUS", "NOT_COMPUTED", "UNSUPPORTED_FORMAT"}
ALLOWED_ASSIGNMENTS = {
    ("VAMOS", "RAW_READS"): ("UNASSIGNED", False),
    ("STRAGLR", "RAW_READS"): ("UNASSIGNED", False),
    ("VAMOS", "ASSEMBLED_CONTIG"): ("DIRECT_SEQUENCE_ASSOCIATION", True),
    ("TANDEM_GENOTYPES", "ASSEMBLED_CONTIG"): ("DIRECT_SEQUENCE_ASSOCIATION", True),
}


@dataclass(frozen=True)
class NormalizationValidationIssue:
    issue_id: str
    severity: str
    code: str
    message: str
    artifact_path: str | None = None
    record_id: str | None = None
    source_file_id: str | None = None
    caller: str | None = None


@dataclass(frozen=True)
class EvidenceSourceReference:
    file_id: str
    caller: str
    caller_version: str | None
    analysis_source: str
    path: str
    sha256: str
    size_bytes: int
    media_type: str
    producer_command_id: str
    associated_sequence_id: str | None


@dataclass(frozen=True)
class CallerEvidenceSummary:
    caller: str
    caller_version: str | None
    analysis_source: str
    run_status: str
    evidence_state: str
    record_count: int
    native_output_count: int
    associated_sequence_ids: tuple[str, ...]
    warnings: tuple[str, ...]
    failure: Mapping[str, Any] | None


class NormalizationError(ValueError):
    """Raised after an integrity failure has been included in a report."""

    def __init__(self, message: str, report: Mapping[str, Any]):
        super().__init__(message)
        self.report = report


def _load(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError("top-level JSON value must be an object")
    return value


def _schema(name: str) -> dict[str, Any]:
    return _load(Path(__file__).resolve().parents[2] / "schemas" / name)


def _issue(issues: list[NormalizationValidationIssue], severity: str, code: str, message: str,
           *, artifact: Path | None = None, record: Mapping[str, Any] | None = None,
           source_id: str | None = None, caller: str | None = None) -> None:
    issues.append(NormalizationValidationIssue(
        issue_id=f"issue-{len(issues)+1:04d}", severity=severity, code=code, message=message,
        artifact_path=str(artifact) if artifact else None,
        record_id=str(record.get("record_id")) if record and record.get("record_id") is not None else None,
        source_file_id=source_id or (str(record.get("source_file_id")) if record and record.get("source_file_id") else None),
        caller=caller or (str(record.get("caller")) if record and record.get("caller") else None)))


def _validate_document(value: Any, schema_name: str, path: Path,
                       issues: list[NormalizationValidationIssue]) -> bool:
    try:
        schema = _schema(schema_name)
        # Caller package schemas refer to the shared record schema by filename.
        # The dependency-free validator intentionally only resolves local refs;
        # validate the envelope here and each record explicitly below.
        items = schema.get("properties", {}).get("records", {}).get("items", {})
        if isinstance(items, dict) and str(items.get("$ref", "")).endswith("normalized-caller-evidence.schema.json"):
            schema = dict(schema)
            schema["properties"] = dict(schema["properties"])
            schema["properties"]["records"] = dict(schema["properties"]["records"])
            schema["properties"]["records"]["items"] = {}
        validate(value, schema)
        return True
    except (SchemaViolation, KeyError, TypeError) as exc:
        _issue(issues, "ERROR", "SCHEMA_VALIDATION_FAILURE", str(exc), artifact=path)
        return False


def _record_key(record: Mapping[str, Any]) -> tuple[str, str, str, str, str, str]:
    # Null sequence IDs sort first and all values are stringified to avoid mixed
    # JSON scalar comparison differences across Python versions.
    return (str(record.get("analysis_source", "")),
            "" if record.get("associated_sequence_id") is None else str(record["associated_sequence_id"]),
            str(record.get("caller", "")), str(record.get("native_record_identifier") or ""),
            str(record.get("native_allele_identifier") or ""), str(record.get("record_id", "")))


def _source_key(source: Mapping[str, Any]) -> tuple[str, str, str, str]:
    return (str(source["caller"]), str(source["analysis_source"]),
            "" if source.get("associated_sequence_id") is None else str(source["associated_sequence_id"]),
            str(source["file_id"]))


def _patient_sequences(metadata: Mapping[str, Any], fasta_path: Path | None,
                       issues: list[NormalizationValidationIssue]) -> dict[str, Mapping[str, Any]]:
    sequences = {str(item["sequence_id"]): item for item in metadata.get("sequences", [])}
    if fasta_path and fasta_path.is_file():
        observed: dict[str, str] = {}
        name: str | None = None
        chunks: list[str] = []
        for line in fasta_path.read_text(encoding="utf-8").splitlines():
            if line.startswith(">"):
                if name is not None: observed[name] = "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else: chunks.append(line.strip())
        if name is not None: observed[name] = "".join(chunks)
        import hashlib
        for sequence_id, item in sequences.items():
            seq = observed.get(sequence_id)
            if seq is None:
                _issue(issues, "ERROR", "PATIENT_SEQUENCE_MISSING", f"{sequence_id} is absent from patient FASTA", artifact=fasta_path)
            elif hashlib.sha256(seq.encode()).hexdigest() != item.get("sequence_sha256"):
                _issue(issues, "ERROR", "PATIENT_SEQUENCE_CHECKSUM_MISMATCH", f"checksum differs for {sequence_id}", artifact=fasta_path)
    return sequences


def normalize_caller_artifacts(*, case_id: str, subject_id: str, sample_id: str, locus_id: str,
                               patient_metadata_path: Path, artifact_groups: Sequence[Mapping[str, Any]],
                               patient_fasta_path: Path | None = None, locus_config_path: Path | None = None,
                               created_utc: str | None = None) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any], dict[str, Any]]:
    """Load and validate adapter products and return package, summary, registry, report.

    Each group supplies ``caller``, ``analysis_source`` and optional ``normalized``,
    ``registry`` and ``run_metadata`` paths. Missing groups are retained as
    NOT_COMPUTED rather than represented as successful empty results.
    """
    del locus_config_path  # Its identity is runner provenance; adapters already used its values.
    issues: list[NormalizationValidationIssue] = []
    metadata = _load(Path(patient_metadata_path))
    _validate_document(metadata, "patient-sequence-metadata.schema.json", Path(patient_metadata_path), issues)
    sequences = _patient_sequences(metadata, Path(patient_fasta_path) if patient_fasta_path else None, issues)
    records: list[dict[str, Any]] = []
    sources: dict[str, dict[str, Any]] = {}
    summaries: list[dict[str, Any]] = []

    supplied = {(str(g["caller"]), str(g["analysis_source"])): g for g in artifact_groups}
    for caller, analysis_source in CALLER_ORDER:
        group = supplied.get((caller, analysis_source), {})
        paths = {key: Path(group[key]) if group.get(key) else None for key in ("normalized", "registry", "run_metadata")}
        present = [p is not None and p.is_file() for p in paths.values()]
        if not all(present):
            if any(present):
                _issue(issues, "ERROR", "INCOMPLETE_CALLER_ARTIFACT_GROUP", f"present {caller}/{analysis_source} artifact group is incomplete", caller=caller)
            else:
                _issue(issues, "WARNING", "MISSING_OPTIONAL_CALLER", f"optional {caller}/{analysis_source} artifacts are absent", caller=caller)
            summaries.append(asdict(CallerEvidenceSummary(caller, None, analysis_source, "NOT_COMPUTED", "NOT_COMPUTED", 0, 0, (), ("Optional caller artifacts are absent.",), None)))
            continue
        try:
            normalized, registry, run = (_load(paths[k]) for k in ("normalized", "registry", "run_metadata"))  # type: ignore[arg-type]
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            _issue(issues, "ERROR", "ARTIFACT_LOAD_FAILURE", str(exc), caller=caller)
            continue
        caller_schema = {"VAMOS": "vamos-normalized-evidence.schema.json", "STRAGLR": "straglr-normalized-evidence.schema.json", "TANDEM_GENOTYPES": "tandem-genotypes-normalized-evidence.schema.json"}[caller]
        valid = _validate_document(normalized, caller_schema, paths["normalized"], issues)  # type: ignore[arg-type]
        valid &= _validate_document(registry, "native-caller-output-registry.schema.json", paths["registry"], issues)  # type: ignore[arg-type]
        # Run schemas differ by adapter but validate their stable outer contract.
        if "runs" in run:
            run_schema = {"VAMOS": "vamos-stage-summary.schema.json", "TANDEM_GENOTYPES": "tandem-genotypes-stage-summary.schema.json"}.get(caller, "straglr-run-metadata.schema.json")
        else:
            run_schema = {"VAMOS": "vamos-run-metadata.schema.json", "STRAGLR": "straglr-run-metadata.schema.json", "TANDEM_GENOTYPES": "tandem-genotypes-run-metadata.schema.json"}[caller]
        valid &= _validate_document(run, run_schema, paths["run_metadata"], issues)  # type: ignore[arg-type]
        if not valid: continue
        if "runs" not in run and (run.get("caller") != caller or run.get("analysis_source") != analysis_source):
            _issue(issues, "ERROR", "CALLER_IDENTITY_MISMATCH", "run metadata caller/source does not match its artifact group", artifact=paths["run_metadata"], caller=caller)
        state = str(normalized.get("evidence_state", "NOT_COMPUTED"))
        if state not in EVIDENCE_STATES:
            _issue(issues, "ERROR", "UNSUPPORTED_EVIDENCE_STATE", f"unsupported evidence state {state}", caller=caller)
        run_status = str(run.get("status", "NOT_COMPUTED"))
        if run_status == "FAILED": state = "COMPUTATION_FAILED"
        caller_version = run.get("caller_version") or next((x.get("caller_version") for x in run.get("runs", []) if x.get("caller_version")), None)
        group_sources: list[dict[str, Any]] = []
        for native in registry.get("outputs", []):
            source = {k: native.get(k) for k in ("file_id", "caller", "caller_version", "analysis_source", "path", "sha256", "size_bytes", "media_type", "producer_command_id")}
            matching_records = [r for r in normalized.get("records", []) if r.get("source_file_id") == native.get("file_id")]
            source["associated_sequence_id"] = native.get("associated_sequence_id", run.get("associated_sequence_id") or next((r.get("associated_sequence_id") for r in matching_records), None))
            file_id = str(source["file_id"])
            existing = sources.get(file_id)
            if existing is not None and existing != source:
                _issue(issues, "ERROR", "SOURCE_REGISTRY_CONFLICT", f"conflicting entries for source file ID {file_id}", artifact=paths["registry"], source_id=file_id, caller=caller)
                continue
            native_path = Path(str(source["path"]))
            if not native_path.is_absolute(): native_path = (paths["registry"].parent / native_path).resolve()  # type: ignore[union-attr]
            if state == "AVAILABLE":
                if not native_path.is_file():
                    _issue(issues, "ERROR", "SOURCE_FILE_MISSING", f"native source does not exist: {native_path}", source_id=file_id, caller=caller)
                elif sha256_file(native_path) != source["sha256"] or native_path.stat().st_size != source["size_bytes"]:
                    _issue(issues, "ERROR", "SOURCE_CHECKSUM_MISMATCH", f"native source identity changed: {native_path}", source_id=file_id, caller=caller)
            sources[file_id] = source
            group_sources.append(source)
        group_records = normalized.get("records", [])
        for record in group_records:
            # Validate shared record independently because the small validator does
            # not dereference cross-file JSON Schema references.
            unified = dict(record)
            unified.setdefault("coordinate_space", "UNKNOWN_COORDINATE_SPACE")
            if not _validate_document(unified, "normalized-caller-evidence.schema.json", paths["normalized"], issues): continue  # type: ignore[arg-type]
            expected = ALLOWED_ASSIGNMENTS.get((record.get("caller"), record.get("analysis_source")))
            if expected is None or record.get("assignment_state") != expected[0] or bool(record.get("associated_sequence_id")) != expected[1]:
                _issue(issues, "ERROR", "INVALID_ASSIGNMENT_COMBINATION", "caller, source, assignment, and sequence association are inconsistent", record=record)
            sequence_id = record.get("associated_sequence_id")
            if sequence_id is not None and sequence_id not in sequences:
                _issue(issues, "ERROR", "UNKNOWN_SEQUENCE_ID", f"unknown patient sequence ID {sequence_id}", record=record)
            for field, expected_identity in (("case_id", case_id), ("subject_id", subject_id), ("sample_id", sample_id), ("locus_id", locus_id)):
                if record.get(field) != expected_identity:
                    _issue(issues, "ERROR", f"{field.upper()}_MISMATCH", f"record {field} does not match package", record=record)
            source = sources.get(str(record.get("source_file_id")))
            if source is None:
                _issue(issues, "ERROR", "UNRESOLVED_SOURCE_FILE", "record source_file_id does not resolve", record=record)
            elif record.get("source_file_sha256") != source.get("sha256"):
                _issue(issues, "ERROR", "SOURCE_CHECKSUM_MISMATCH", "record checksum differs from source registry", record=record)
            records.append(unified)
        warnings = tuple(str(x) for x in (*run.get("warnings", []), *normalized.get("normalization_warnings", [])))
        summaries.append(asdict(CallerEvidenceSummary(caller, caller_version, analysis_source, run_status, state,
            len(group_records), len(group_sources), tuple(sorted({str(r["associated_sequence_id"]) for r in group_records if r.get("associated_sequence_id")})), warnings, run.get("failure"))))

    seen_records: set[str] = set(); seen_native: set[tuple[Any, ...]] = set()
    for record in records:
        if record["record_id"] in seen_records: _issue(issues, "ERROR", "DUPLICATE_RECORD_ID", "duplicate unified record ID", record=record)
        seen_records.add(record["record_id"])
        native_key = (record["caller"], record["analysis_source"], record.get("native_record_identifier"), record.get("native_allele_identifier"))
        if native_key in seen_native: _issue(issues, "ERROR", "DUPLICATE_NATIVE_IDENTIFIER", "duplicate caller/native record/native allele identity", record=record)
        seen_native.add(native_key)
    records.sort(key=_record_key)
    source_list = sorted(sources.values(), key=_source_key)
    summaries.sort(key=lambda x: CALLER_ORDER.index((x["caller"], x["analysis_source"])))
    groups = []
    for sequence_id in (None, *sorted(sequences)):
        applicable = [s for s in summaries if (sequence_id is None) == (s["analysis_source"] == "RAW_READS")]
        counts = {f'{r["caller"]}:{r["analysis_source"]}': sum(1 for item in records if item["caller"] == r["caller"] and item["analysis_source"] == r["analysis_source"] and item.get("associated_sequence_id") == sequence_id) for r in applicable}
        groups.append({"sequence_id": sequence_id, "available_callers": [f'{s["caller"]}:{s["analysis_source"]}' for s in applicable if s["evidence_state"] == "AVAILABLE"],
                       "missing_callers": [f'{s["caller"]}:{s["analysis_source"]}' for s in applicable if s["evidence_state"] in {"NOT_COMPUTED", "INPUT_MISSING", "NOT_APPLICABLE"}],
                       "unsupported_callers": [f'{s["caller"]}:{s["analysis_source"]}' for s in applicable if s["evidence_state"] == "UNSUPPORTED_FORMAT"],
                       "failed_callers": [f'{s["caller"]}:{s["analysis_source"]}' for s in applicable if s["evidence_state"] == "COMPUTATION_FAILED"], "record_counts": counts})
    created = created_utc or utc_now()
    package_core = {"record_schema_version": SCHEMA_VERSION, "case_id": case_id, "subject_id": subject_id,
                    "sample_id": sample_id, "locus_id": locus_id, "records": records,
                    "source_registry": source_list, "caller_summaries": summaries,
                    "normalization_warnings": sorted({w for s in summaries for w in s["warnings"]})}
    package = {**package_core, "package_id": f"normalized-{canonical_digest(package_core)[:24]}", "created_utc": created}
    summary = {"record_schema_version": SCHEMA_VERSION, "package_id": package["package_id"], "case_id": case_id,
               "caller_summaries": summaries, "comparison_availability": groups,
               "note": "Inventory only; shared fields do not establish equivalent measurements."}
    registry_document = {"record_schema_version": SCHEMA_VERSION, "sources": source_list}
    report = {"record_schema_version": SCHEMA_VERSION, "package_id": package["package_id"], "valid": not any(i.severity == "ERROR" for i in issues),
              "error_count": sum(i.severity == "ERROR" for i in issues), "warning_count": sum(i.severity == "WARNING" for i in issues),
              "issues": [asdict(i) for i in issues]}
    if not report["valid"]: raise NormalizationError("unified normalization validation failed", report)
    return package, summary, registry_document, report


def write_normalization_outputs(output_directory: Path, documents: tuple[Mapping[str, Any], ...], *, overwrite: bool = False) -> dict[str, Path]:
    """Atomically publish the four final contract documents."""
    output_directory = Path(output_directory)
    if output_directory.exists() and not overwrite: raise FileExistsError(f"output already exists: {output_directory}")
    output_directory.parent.mkdir(parents=True, exist_ok=True)
    work = output_directory.parent / f".work-{uuid.uuid4()}"
    work.mkdir()
    names = ("normalized-evidence.json", "evidence-summary.json", "source-registry.json", "validation-report.json")
    try:
        for name, document in zip(names, documents): atomic_write_json(work / name, document)
        if output_directory.exists(): shutil.rmtree(output_directory)
        os.replace(work, output_directory)
    except Exception:
        shutil.rmtree(work, ignore_errors=True)
        raise
    return {name: output_directory / name for name in names}


def run_normalization_stage(config: Mapping[str, Any], root: Path, *, overwrite: bool = False) -> dict[str, Path]:
    run = config["run"]
    specs = (
        ("VAMOS", "RAW_READS", "04_vamos_read", "vamos-read.normalized.json", "vamos-read.outputs.json", "vamos-read.run.json"),
        ("VAMOS", "ASSEMBLED_CONTIG", "05_vamos_contig", "stage-normalized.json", "stage-outputs.json", "stage-summary.json"),
        ("STRAGLR", "RAW_READS", "06_straglr", "straglr.normalized.json", "straglr.outputs.json", "straglr.run.json"),
        ("TANDEM_GENOTYPES", "ASSEMBLED_CONTIG", "08_tandem_genotypes", "stage-normalized.json", "stage-outputs.json", "stage-summary.json"),
    )
    groups = [{"caller": c, "analysis_source": s, "normalized": root / d / n,
               "registry": root / d / o, "run_metadata": root / d / m} for c, s, d, n, o, m in specs]
    try:
        documents = normalize_caller_artifacts(case_id=run["case_id"], subject_id=run["subject_id"], sample_id=run["sample_id"],
            locus_id=run["locus_id"], patient_metadata_path=root / "03_assembly_alignment" / "patient-sequences.metadata.json",
            patient_fasta_path=root / "03_assembly_alignment" / "patient-sequences.fasta", locus_config_path=Path(str(config["locus_config"])), artifact_groups=groups)
    except NormalizationError as exc:
        # A report is diagnostic, not a biological package. Publish only that
        # report on failure; no normalized package, summary, or registry can be
        # mistaken for a successful partial result.
        failed = root / "09_normalized_evidence"
        if failed.exists() and not overwrite:
            raise
        if failed.exists(): shutil.rmtree(failed)
        failed.mkdir(parents=True)
        atomic_write_json(failed / "validation-report.json", exc.report)
        raise
    return write_normalization_outputs(root / "09_normalized_evidence", documents, overwrite=overwrite)
