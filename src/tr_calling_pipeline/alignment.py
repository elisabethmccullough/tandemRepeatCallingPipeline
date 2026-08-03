"""SAM alignment classification and stable mapping/orientation contracts."""
from __future__ import annotations
from dataclasses import asdict, dataclass
import re
from typing import Iterable

from .sequence_metadata import AlignmentStrand, DisplayOrientation, MappingStatus

_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def cigar_lengths(cigar: str) -> tuple[int, int]:
    if cigar == "*": return 0, 0
    fields = _CIGAR.findall(cigar)
    if not fields or "".join(n + op for n, op in fields) != cigar:
        raise ValueError(f"invalid CIGAR: {cigar}")
    reference = sum(int(n) for n, op in fields if op in "MDN=X")
    mapped = sum(int(n) for n, op in fields if op in "M=X")
    return reference, mapped


@dataclass(frozen=True)
class SamAlignment:
    query: str; flag: int; contig: str | None; start: int | None; mapq: int | None; cigar: str | None
    @property
    def unmapped(self): return bool(self.flag & 4)
    @property
    def secondary(self): return bool(self.flag & 256)
    @property
    def supplementary(self): return bool(self.flag & 2048)


@dataclass(frozen=True)
class AssemblyAlignmentMapping:
    record_schema_version: str
    sequence_id: str
    source_fasta_record_id: str
    sequence_sha256: str
    reference_contig: str | None
    alignment_start: int | None
    alignment_end: int | None
    coordinate_convention: str
    reference_alignment_strand: str
    mapping_quality: int | None
    cigar: str | None
    primary_alignment_count: int
    secondary_alignment_count: int
    supplementary_alignment_count: int
    mapped_base_count: int | None
    source_sequence_length: int
    alignment_span: int | None
    original_orientation: str
    display_orientation: str
    reverse_complement_required: bool | None
    mapping_status: str
    warnings: tuple[str, ...]
    def to_dict(self):
        value = asdict(self); value["warnings"] = list(self.warnings); return value


def parse_sam(lines: Iterable[str]) -> tuple[SamAlignment, ...]:
    result=[]
    for number, line in enumerate(lines, 1):
        if not line.strip() or line.startswith("@"): continue
        fields=line.rstrip("\n").split("\t")
        if len(fields) < 11: raise ValueError(f"invalid SAM record at line {number}")
        try: flag, pos, mapq = int(fields[1]), int(fields[3]), int(fields[4])
        except ValueError as exc: raise ValueError(f"invalid SAM numeric field at line {number}") from exc
        unmapped=bool(flag & 4)
        result.append(SamAlignment(fields[0], flag, None if unmapped or fields[2]=="*" else fields[2],
            None if unmapped or pos == 0 else pos-1, None if unmapped else mapq, None if fields[5]=="*" else fields[5]))
    return tuple(result)


def classify_mapping(sequence_id: str, source_record_id: str, digest: str, length: int, alignments: Iterable[SamAlignment]) -> AssemblyAlignmentMapping:
    items=[a for a in alignments if a.query == sequence_id]
    mapped=[a for a in items if not a.unmapped]
    primary=[a for a in mapped if not a.secondary and not a.supplementary]
    secondary=[a for a in mapped if a.secondary and not a.supplementary]
    supplementary=[a for a in mapped if a.supplementary]
    chosen=primary[0] if len(primary)==1 else None
    if len(primary)>1: status=MappingStatus.MULTIPLE_PRIMARY.value
    elif not mapped: status=MappingStatus.UNMAPPED.value
    elif not primary and secondary and not supplementary: status=MappingStatus.SECONDARY_ONLY.value
    elif not primary and supplementary and not secondary: status=MappingStatus.SUPPLEMENTARY_ONLY.value
    elif not primary: status=MappingStatus.AMBIGUOUS.value
    else: status=MappingStatus.UNIQUE_PRIMARY.value
    span=mapped_bases=None
    if chosen and chosen.cigar: span,mapped_bases=cigar_lengths(chosen.cigar)
    strand=(AlignmentStrand.REVERSE.value if chosen and chosen.flag & 16 else AlignmentStrand.FORWARD.value) if status==MappingStatus.UNIQUE_PRIMARY.value else AlignmentStrand.UNKNOWN.value
    resolved=status==MappingStatus.UNIQUE_PRIMARY.value
    warnings=() if resolved else (f"mapping status is {status}",)
    return AssemblyAlignmentMapping("1.0",sequence_id,source_record_id,digest,chosen.contig if chosen else None,
        chosen.start if chosen else None,(chosen.start+span) if chosen and chosen.start is not None and span is not None else None,
        "zero_based_half_open",strand,chosen.mapq if chosen else None,chosen.cigar if chosen else None,len(primary),len(secondary),len(supplementary),mapped_bases,length,span,
        "AS_STORED",DisplayOrientation.REFERENCE_FORWARD.value if resolved else DisplayOrientation.UNRESOLVED.value,
        (strand==AlignmentStrand.REVERSE.value) if resolved else None,status,warnings)


def align_assembly(config: dict, output_directory, minimap2_tool, samtools_tool, *, overwrite=False) -> dict:
    """Run a deliberately unpiped alignment flow, preserving one record per command."""
    import json
    from pathlib import Path
    from .execution import CommandSpec, InputDeclaration, OutputDeclaration, execute
    from .fasta import write_fasta
    from .provenance import atomic_write_json, sha256_file
    from .sequence_metadata import select_patient_sequences, with_mapping

    stage=Path(output_directory).resolve(); stage.mkdir(parents=True,exist_ok=True)
    if not minimap2_tool.resolved_executable or not samtools_tool.resolved_executable:
        raise ValueError("required minimap2 and samtools executables must be available")
    reference=Path(config["inputs"]["reference_fasta"]); assembly=Path(config["inputs"]["assembly_fasta"])
    selected,metadata=select_patient_sequences(assembly,config["inputs"]["assembly_records"])
    package_fasta=stage/"patient-sequences.fasta"
    if package_fasta.exists() and not overwrite: raise ValueError(f"output already exists: {package_fasta}")
    write_fasta(package_fasta,selected)
    settings={"minimap2_preset":"asm20","secondary":True,"threads":config["execution"]["threads"]}
    settings.update(config.get("assembly_alignment",{}))
    records=stage/"execution-records"; logs=stage/"logs"; records.mkdir(exist_ok=True)
    def run(tool,cid,argv,inputs,outputs=()):
        return execute(CommandSpec(cid,"02_align_assembly",tool.tool_id.value,tuple(argv),str(stage),
            declared_inputs=tuple(InputDeclaration(str(i),str(p)) for i,p in enumerate(inputs)),
            declared_outputs=tuple(OutputDeclaration(str(i),str(p)) for i,p in enumerate(outputs)),overwrite=overwrite),tool,records/f"{cid}.json",logs)
    secondary="yes" if settings["secondary"] else "no"
    mm=run(minimap2_tool,"minimap2",(str(minimap2_tool.resolved_executable),"-a","-x",str(settings["minimap2_preset"]),"-t",str(settings["threads"]),"--secondary="+secondary,str(reference),str(package_fasta)),(reference,package_fasta))
    sam=stage/"assembly.aligned.sam"; sam.write_bytes(Path(mm["stdout_log"]).read_bytes())
    unsorted=stage/"assembly.aligned.unsorted.bam"; sorted_bam=stage/"assembly.aligned.sorted.bam"; index=Path(str(sorted_bam)+".bai")
    run(samtools_tool,"samtools-view",(str(samtools_tool.resolved_executable),"view","-b","-o",str(unsorted),str(sam)),(sam,), (unsorted,))
    run(samtools_tool,"samtools-sort",(str(samtools_tool.resolved_executable),"sort","-o",str(sorted_bam),str(unsorted)),(unsorted,), (sorted_bam,))
    run(samtools_tool,"samtools-index",(str(samtools_tool.resolved_executable),"index","-o",str(index),str(sorted_bam)),(sorted_bam,), (index,))
    flag=run(samtools_tool,"flagstat",(str(samtools_tool.resolved_executable),"flagstat",str(sorted_bam)),(sorted_bam,index))
    idx=run(samtools_tool,"idxstats",(str(samtools_tool.resolved_executable),"idxstats",str(sorted_bam)),(sorted_bam,index))
    (stage/"assembly_alignment.flagstat.txt").write_bytes(Path(flag["stdout_log"]).read_bytes())
    (stage/"assembly_alignment.idxstats.txt").write_bytes(Path(idx["stdout_log"]).read_bytes())
    alignments=parse_sam(sam.read_text().splitlines())
    mappings=[classify_mapping(m.sequence_id,m.source_fasta_record_id,m.sequence_sha256,m.sequence_length,alignments) for m in metadata]
    mapping_doc={"schema_version":"1.0","coordinate_convention":"zero_based_half_open","mappings":[m.to_dict() for m in mappings]}
    atomic_write_json(stage/"assembly_record_mappings.json",mapping_doc)
    updated=[]
    for meta,mapping in zip(metadata,mappings):
        coordinates=None if mapping.alignment_start is None else {"reference_contig":mapping.reference_contig,"start":mapping.alignment_start,"end":mapping.alignment_end,"coordinate_convention":"zero_based_half_open"}
        updated.append(with_mapping(meta,strand=mapping.reference_alignment_strand,status=mapping.mapping_status,coordinates=coordinates,warnings=mapping.warnings).to_dict())
    atomic_write_json(stage/"patient-sequences.metadata.json",{"schema_version":"1.0","patient_fasta_path":str(package_fasta),"patient_fasta_sha256":sha256_file(package_fasta),"sequences":updated})
    summary={"schema_version":"1.0","alignment_configuration":settings,"development_default_warning":"asm20 is a development default pending biological validation.",
        "source_assembly_fasta_sha256":sha256_file(assembly),"package_patient_fasta_sha256":sha256_file(package_fasta),"aligned_bam_sha256":sha256_file(sorted_bam),"aligned_bam_index_sha256":sha256_file(index),
        "sequence_count":len(metadata),"mapping_status_counts":{status.value:sum(m.mapping_status==status.value for m in mappings) for status in MappingStatus if status is not MappingStatus.NOT_ALIGNED},"warnings":[]}
    atomic_write_json(stage/"assembly_alignment.summary.json",summary)
    return summary
