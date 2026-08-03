"""Offline synthetic demonstrations with explicit fake-tool boundaries."""
from __future__ import annotations

import hashlib
from importlib.resources import files
import json
import os
from pathlib import Path
import shutil
import stat
import sys

import yaml

from .case_package import build_case_package
from .case_package_validation import validate_case_package
from .normalization import run_normalization_stage
from .provenance import sha256_file
from .runner import run

EXPECTED_GROUPS = (("VAMOS", "RAW_READS"), ("VAMOS", "ASSEMBLED_CONTIG"),
                   ("STRAGLR", "RAW_READS"), ("TANDEM_GENOTYPES", "ASSEMBLED_CONTIG"))


def run_package_portability_demo(output: str | Path) -> dict[str, object]:
    """Focused package-only check retained separately from the runner demo."""
    output = Path(output).resolve()
    if output.exists(): raise FileExistsError(f"demo output already exists: {output}")
    source = output / "temporary package-only run"
    patient = source / "03_assembly_alignment"; patient.mkdir(parents=True)
    records = [("synthetic-sequence-1", "ACGTACGT"), ("synthetic-sequence-2", "TGCATGCA")]
    (patient / "patient-sequences.fasta").write_text("".join(f">{name}\n{seq}\n" for name, seq in records))
    sequences=[]
    for index,(name,sequence) in enumerate(records,1):
        sequences.append({"record_schema_version":"1.0","sequence_id":name,"sequence_role":"PATIENT_HAPLOTYPE",
            "display_label":f"Synthetic sequence {index}","source_fasta_record_id":f"SYNTHETIC_{index}",
            "sequence_sha256":hashlib.sha256(sequence.encode()).hexdigest(),"sequence_length":len(sequence),
            "original_orientation":"AS_STORED","reference_alignment_strand":"NOT_ALIGNED","display_orientation":"UNRESOLVED",
            "reverse_complement_required":None,"source_coordinates":None,"mapping_status":"NOT_ALIGNED","warnings":["SYNTHETIC_NON_CLINICAL"]})
    (patient/"patient-sequences.metadata.json").write_text(json.dumps({"schema_version":"1.0","sequences":sequences}))
    (patient/"assembly_record_mappings.json").write_text(json.dumps({"schema_version":"1.0","coordinate_convention":"zero_based_half_open","mappings":[]}))
    controls=output/"package-only controls"; controls.mkdir(parents=True)
    locus=controls/"locus.yaml"; locus.write_text('schema_version: "1.0"\n')
    config_path=controls/"run.yaml"; config_path.write_text('schema_version: "1.0"\n')
    identity={"case_id":"SYNTHETIC-PORTABILITY","subject_id":"SYNTHETIC-SUBJECT","sample_id":"SYNTHETIC-SAMPLE","locus_id":"SYNTHETIC-LOCUS"}
    run_normalization_stage({"run":identity,"locus_config":str(locus)},source)
    built=output/"built package"
    build_case_package({"run":identity,"locus_config":str(locus),"_config_path":str(config_path),"case_package":{
        "package_root":str(built),"include_native_outputs":True,"include_prepared_mini_bam":False,
        "include_alignment_artifacts":False,"include_command_records":False}},source)
    before=validate_case_package(built); moved=output/"portable package-only case"; shutil.move(built,moved)
    after=validate_case_package(moved); shutil.rmtree(source); shutil.rmtree(controls)
    removed=validate_case_package(moved,write_report=True)
    return {"package":str(moved),"valid_before_move":before["valid"],"valid_after_move":after["valid"],
            "valid_after_source_removal":removed["valid"],"scope":"PACKAGE_PORTABILITY_ONLY"}


def _fallback_tool(directory: Path, tool: str) -> str:
    """Create a POSIX source-checkout launcher; installed wheels use console scripts."""
    if os.name == "nt":
        raise FileNotFoundError(f"installed fake-{tool} console script is required on Windows")
    path=directory/f"fake-{tool}"
    package_root=Path(__file__).resolve().parents[1]
    path.write_text(f"#!{sys.executable}\nimport sys\nsys.path.insert(0,{str(package_root)!r})\n"
        f"from tr_calling_pipeline.fake_tools import {tool.replace('-', '_')}\nsys.exit({tool.replace('-', '_')}())\n")
    path.chmod(path.stat().st_mode|stat.S_IXUSR)
    return str(path)


def _prepare_workspace(output: Path) -> tuple[Path, Path, dict[str, str]]:
    workspace=output/"synthetic runner workspace"; fixtures=workspace/"fixtures"; fixtures.mkdir(parents=True)
    resource=files("tr_calling_pipeline").joinpath("demo_fixtures")
    for item in resource.iterdir():
        if item.is_file(): shutil.copyfile(str(item),fixtures/item.name)
    launchers=workspace/"fake tool launchers"; launchers.mkdir()
    tools={}
    for name in ("samtools","minimap2","vamos","straglr","lastdb","lastal","tandem-genotypes"):
        tools[name]=shutil.which(f"fake-{name}") or _fallback_tool(launchers,name)
    package=output/"initial synthetic case package"
    config={"schema_version":"1.0","run":{"case_id":"SYNTHETIC-DEMO","subject_id":"SYNTHETIC-SUBJECT",
        "sample_id":"SYNTHETIC-SAMPLE","locus_id":"SYNTHETIC-LOCUS","output_root":str(workspace/"run output"),"overwrite":False},
        "inputs":{"assembly_fasta":str(fixtures/"assembly.fasta"),"assembly_records":[
            {"record_id":"SYNTHETIC_A","sequence_id":"synthetic-sequence-1","display_label":"Synthetic sequence 1","sequence_role":"PATIENT_HAPLOTYPE"},
            {"record_id":"SYNTHETIC_B","sequence_id":"synthetic-sequence-2","display_label":"Synthetic sequence 2","sequence_role":"PATIENT_HAPLOTYPE"}],
            "mini_bam":str(fixtures/"mini.bam"),"mini_bam_index":str(fixtures/"mini.bam.bai"),
            "reference_fasta":str(fixtures/"reference.fasta"),"reference_fasta_index":str(fixtures/"reference.fasta.fai"),
            "reference_scope":"LOCAL_LOCUS_REFERENCE"},"locus_config":str(fixtures/"locus.yaml"),
        "case_package":{"enabled":True,"package_root":str(package),"include_prepared_mini_bam":True,
            "include_alignment_artifacts":True,"include_native_outputs":True,"include_command_records":True,
            "comparison_panel":{"panel_id":"SYNTHETIC_PANEL","panel_version":"development","required":False}},
        "execution":{"threads":1,"dry_run":False},"tools":{
            "samtools":{"executable":tools["samtools"],"required":True,"execution_mode":"NATIVE"},
            "minimap2":{"executable":tools["minimap2"],"required":True,"execution_mode":"NATIVE"},
            "vamos":{"executable":tools["vamos"],"required":True,"execution_mode":"NATIVE","allow_provisional_adapter":True,"additional_arguments":[]},
            "straglr":{"executable":tools["straglr"],"required":True,"execution_mode":"NATIVE"},
            "lastdb":{"executable":tools["lastdb"],"required":True,"execution_mode":"NATIVE"},
            "lastal":{"executable":tools["lastal"],"required":True,"execution_mode":"NATIVE"},
            "tandem_genotypes":{"executable":tools["tandem-genotypes"],"required":True,"execution_mode":"NATIVE"}}}
    config_path=workspace/"synthetic run config.yaml"; config_path.write_text(yaml.safe_dump(config,sort_keys=False))
    return config_path,package,tools


def run_demo(output: str | Path) -> dict[str, object]:
    """Execute all eleven implemented stages through the canonical runner."""
    output=Path(output).resolve()
    if output.exists(): raise FileExistsError(f"demo output already exists: {output}")
    output.mkdir(parents=True)
    config_path,package,tools=_prepare_workspace(output)
    fixture_dir=config_path.parent/"fixtures"
    before={p.name:sha256_file(p) for p in fixture_dir.iterdir() if p.is_file()}
    run_root=Path(run(config_path))
    stage_records={p.stem:json.loads(p.read_text()) for p in sorted((run_root/"00_manifest/stages").glob("*.json")) if ".invalidated." not in p.name}
    if set(stage_records) != {f"{index:02d}_{name}" for index,name in enumerate(("validate_inputs","prepare_bam","align_assembly","run_vamos_read","run_vamos_contig","run_straglr","prepare_tandem_genotypes","run_tandem_genotypes","normalize_outputs","build_case_package","validate_case_package"))}:
        raise RuntimeError("the canonical eleven-stage record set was not produced")
    if any(record["status"]!="SUCCEEDED" for record in stage_records.values()): raise RuntimeError("not every pipeline stage succeeded")
    evidence=json.loads((run_root/"09_normalized_evidence/normalized-evidence.json").read_text())
    counts={f"{caller} / {source}":sum(r["caller"]==caller and r["analysis_source"]==source for r in evidence["records"])
            for caller,source in EXPECTED_GROUPS}
    if any(value<=0 for value in counts.values()): raise RuntimeError(f"missing synthetic evidence group: {counts}")
    tool_records=[tool for record in stage_records.values() for tool in record["tool_identities"] if tool["tool_id"] not in {"PIPELINE","PYTHON"}]
    if any(tool["executable_kind"]!="FAKE_TEST_TOOL" or tool["verification_level"]!="SYNTHETIC_INTEGRATION_TESTED" for tool in tool_records):
        raise RuntimeError("synthetic executable provenance is not explicitly fake and synthetic-tested")
    after={p.name:sha256_file(p) for p in fixture_dir.iterdir() if p.is_file()}
    if before != after: raise RuntimeError("a synthetic source fixture was modified")
    before_move=validate_case_package(package)
    # An unchanged second canonical run must conservatively resume every stage.
    run(config_path)
    skip_count=len(list((run_root/"00_manifest/stages").glob("*.skip.json")))
    moved=output/"portable synthetic case package"; shutil.move(package,moved)
    after_move=validate_case_package(moved)
    shutil.rmtree(config_path.parent)
    after_removal=validate_case_package(moved,write_report=True)
    return {"synthetic":True,"non_clinical":True,"real_tools_executed":False,
        "verification_level":"SYNTHETIC_INTEGRATION_TESTED","tool_notice":"Only bundled FAKE_TEST_TOOL executables ran; no real caller or laboratory workflow was verified.",
        "package":str(moved),"caller_records_by_source":counts,"stage_count":len(stage_records),"resumed_stage_count":skip_count,
        "valid_before_move":before_move["valid"],"valid_after_move":after_move["valid"],"valid_after_source_removal":after_removal["valid"]}
