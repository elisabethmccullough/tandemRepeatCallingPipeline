"""Samtools-backed mini-BAM validation and immutable-source preparation."""
from __future__ import annotations
from dataclasses import asdict, dataclass
import errno, os, re, shutil
from pathlib import Path
from typing import Any

from .execution import CommandSpec, InputDeclaration, OutputDeclaration, execute
from .errors import OutputConflictError
from .provenance import atomic_write_json, sha256_file
from .tools import Tool


@dataclass(frozen=True)
class BamCountSummary:
    record_schema_version: str = "1.0"
    total_reads: int | None = None
    primary_reads: int | None = None
    secondary_reads: int | None = None
    supplementary_reads: int | None = None
    mapped_reads: int | None = None
    unmapped_reads: int | None = None
    duplicate_reads: int | None = None
    paired_reads: int | None = None
    properly_paired_reads: int | None = None


@dataclass(frozen=True)
class BamValidationResult:
    schema_version: str
    valid: bool
    bam_path: str
    index_path: str
    sort_order: str
    reference_contigs: tuple[str, ...]
    counts: BamCountSummary
    warnings: tuple[str, ...] = ()


def parse_header(text: str) -> tuple[str, tuple[str, ...]]:
    sort_order="unknown"; contigs=[]
    for line in text.splitlines():
        fields=line.split("\t")
        if fields and fields[0]=="@HD":
            for field in fields[1:]:
                if field.startswith("SO:"): sort_order=field[3:]
        elif fields and fields[0]=="@SQ":
            for field in fields[1:]:
                if field.startswith("SN:"): contigs.append(field[3:]); break
    return sort_order,tuple(contigs)


def parse_flagstat(text: str) -> BamCountSummary:
    values: dict[str,int] = {}
    patterns={"total_reads":r"^(\d+) \+ \d+ in total", "secondary_reads":r"^(\d+) \+ \d+ secondary",
        "supplementary_reads":r"^(\d+) \+ \d+ supplementary", "duplicate_reads":r"^(\d+) \+ \d+ duplicates",
        "mapped_reads":r"^(\d+) \+ \d+ mapped", "paired_reads":r"^(\d+) \+ \d+ paired in sequencing",
        "properly_paired_reads":r"^(\d+) \+ \d+ properly paired"}
    for line in text.splitlines():
        for name,pattern in patterns.items():
            match=re.match(pattern,line)
            if match: values[name]=int(match.group(1))
    total=values.get("total_reads"); secondary=values.get("secondary_reads"); supplementary=values.get("supplementary_reads")
    primary=(total-secondary-supplementary) if None not in (total,secondary,supplementary) else None
    mapped=values.get("mapped_reads")
    return BamCountSummary(total_reads=total,primary_reads=primary,secondary_reads=secondary,
        supplementary_reads=supplementary,mapped_reads=mapped,unmapped_reads=total-mapped if total is not None and mapped is not None else None,
        duplicate_reads=values.get("duplicate_reads"),paired_reads=values.get("paired_reads"),properly_paired_reads=values.get("properly_paired_reads"))


def _run(tool: Tool, stage: Path, command_id: str, args: tuple[str,...], bam: Path,
         index: Path | None = None, *, output: Path | None = None, overwrite: bool = False):
    record=stage/"execution-records"/f"{command_id}.json"; record.parent.mkdir(parents=True,exist_ok=True)
    declared=() if output is None else (OutputDeclaration(command_id,str(output),role="BAM"),)
    inputs = [InputDeclaration("bam", str(bam), role="BAM")]
    if index is not None:
        inputs.append(InputDeclaration("index", str(index), role="BAM_INDEX"))
    spec=CommandSpec(command_id,"01_prepare_bam",tool.tool_id.value,(str(tool.resolved_executable),*args),str(stage.resolve()),
        declared_inputs=tuple(inputs),declared_outputs=declared,overwrite=overwrite)
    return execute(spec,tool,record,stage/"logs")


_LINK_FALLBACK_ERRNOS = frozenset({errno.EXDEV, errno.EPERM, errno.EOPNOTSUPP, errno.ENOSYS})


def _prepare_destination(destination: Path, *, overwrite: bool) -> None:
    if not destination.exists() and not destination.is_symlink():
        return
    if not overwrite:
        raise OutputConflictError(
            f"output already exists and overwrite is disabled: {destination}",
            output_path=str(destination),
        )
    if destination.is_dir() and not destination.is_symlink():
        raise OutputConflictError(f"output path is a directory and cannot be replaced: {destination}")
    destination.unlink()


def _prepare_destinations(destinations: tuple[Path, ...], *, overwrite: bool) -> None:
    """Validate every final output before removing any of them."""
    conflicts = [path for path in destinations if path.exists() or path.is_symlink()]
    if conflicts and not overwrite:
        raise OutputConflictError(
            "output already exists and overwrite is disabled: " + ", ".join(map(str, conflicts)),
            output_paths=[str(path) for path in conflicts],
        )
    directories = [path for path in conflicts if path.is_dir() and not path.is_symlink()]
    if directories:
        raise OutputConflictError(
            "output path is a directory and cannot be replaced: " + ", ".join(map(str, directories))
        )
    for path in conflicts:
        path.unlink()


def _link_or_copy(source: Path, destination: Path, *, overwrite: bool = False) -> str:
    """Hard-link when supported, with explicit conflicts and narrow copy fallback."""
    _prepare_destination(destination, overwrite=overwrite)
    try:
        os.link(source, destination)
        return "HARD_LINK"
    except OSError as exc:
        if exc.errno not in _LINK_FALLBACK_ERRNOS:
            raise
    # copy2 is reached only for a known cross-filesystem/capability failure.
    shutil.copy2(source, destination)
    return "COPY"


def _stage_configured_pair(source_bam: Path, source_index: Path, stage: Path) -> tuple[Path, Path, str, str]:
    """Isolate the exact configured pair under names samtools will associate."""
    validation_dir = stage / ".configured-index-validation"
    if validation_dir.exists():
        shutil.rmtree(validation_dir)
    validation_dir.mkdir(parents=True)
    staged_bam = validation_dir / "configured-source.bam"
    staged_index = validation_dir / "configured-source.bam.bai"
    bam_strategy = _link_or_copy(source_bam, staged_bam)
    index_strategy = _link_or_copy(source_index, staged_index)
    return staged_bam, staged_index, bam_strategy, index_strategy


def prepare_bam(source_bam: str|Path, source_index: str|Path, output_directory: str|Path, tool: Tool, *, overwrite=False) -> dict[str,Any]:
    source_bam,source_index=Path(source_bam).resolve(),Path(source_index).resolve(); stage=Path(output_directory).resolve()
    for source in (source_bam,source_index):
        if not source.is_file() or source.stat().st_size==0: raise ValueError(f"BAM input is missing or empty: {source}")
    if not tool.resolved_executable: raise ValueError("required samtools executable is unavailable")
    stage.mkdir(parents=True,exist_ok=True)
    prepared=stage/"prepared.mini.bam"; prepared_index=stage/"prepared.mini.bam.bai"
    # Fail before validation commands can mutate stage state or copy2 can replace data.
    _prepare_destinations((prepared, prepared_index), overwrite=overwrite)
    staged_bam, staged_index, staged_bam_strategy, staged_index_strategy = _stage_configured_pair(
        source_bam, source_index, stage
    )
    _run(tool,stage,"quickcheck-source",("quickcheck",str(staged_bam)),staged_bam)
    _run(tool,stage,"header-source",("view","-H",str(staged_bam)),staged_bam)
    _run(tool,stage,"flagstat-source",("flagstat",str(staged_bam)),staged_bam)
    _run(tool,stage,"idxstats-source",("idxstats",str(staged_bam)),staged_bam,staged_index)
    logs=stage/"logs"/"01_prepare_bam"
    header=(logs/"header-source.stdout.log").read_text(); flagstat=(logs/"flagstat-source.stdout.log").read_text(); idxstats=(logs/"idxstats-source.stdout.log").read_text()
    sort_order,contigs=parse_header(header)
    (stage/"input_bam.header.sam").write_text(header); (stage/"input_bam.flagstat.txt").write_text(flagstat); (stage/"input_bam.idxstats.txt").write_text(idxstats)
    if sort_order=="coordinate":
        strategy=_link_or_copy(source_bam,prepared); index_strategy=_link_or_copy(source_index,prepared_index)
    else:
        temporary=stage/".prepared.mini.bam.tmp"; temporary.unlink(missing_ok=True)
        _run(tool,stage,"sort",("sort","-o",str(temporary),str(staged_bam)),staged_bam,output=temporary,overwrite=True)
        _run(tool,stage,"quickcheck-sorted",("quickcheck",str(temporary)),temporary)
        os.replace(temporary,prepared)
        _run(tool,stage,"index-prepared",("index","-o",str(prepared_index),str(prepared)),prepared,output=prepared_index,overwrite=overwrite)
        strategy="SORTED"; index_strategy="GENERATED"
    # Prove the prepared index can service the prepared BAM.
    _run(tool,stage,"quickcheck-prepared",("quickcheck",str(prepared)),prepared)
    _run(tool,stage,"idxstats-prepared",("idxstats",str(prepared)),prepared,prepared_index)
    validation=BamValidationResult("1.0",True,str(source_bam),str(source_index),sort_order,contigs,parse_flagstat(flagstat))
    atomic_write_json(stage/"input_bam.validation.json",{**asdict(validation),"reference_contigs":list(contigs),"warnings":[]})
    metadata={"schema_version":"1.0","source_bam_sha256":sha256_file(source_bam),"source_index_sha256":sha256_file(source_index),
        "prepared_bam_sha256":sha256_file(prepared),"prepared_index_sha256":sha256_file(prepared_index),"preparation_strategy":strategy,
        "index_preparation_strategy":index_strategy,"source_sort_order":sort_order,"prepared_sort_order":"coordinate","samtools_version":tool.detected_version,
        "configured_source_bam":str(source_bam),"configured_source_index":str(source_index),
        "staged_validation_bam":str(staged_bam),"staged_validation_index":str(staged_index),
        "staged_validation_bam_sha256":sha256_file(staged_bam),"staged_validation_index_sha256":sha256_file(staged_index),
        "staged_validation_bam_strategy":staged_bam_strategy,"staged_validation_index_strategy":staged_index_strategy,
        "command_record_paths":[str(path) for path in sorted((stage/"execution-records").glob("*.json"))]}
    atomic_write_json(stage/"prepared_bam.metadata.json",metadata)
    return metadata
