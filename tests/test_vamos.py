import json
from pathlib import Path
import shutil
import pytest
import yaml
from tr_calling_pipeline.callers.vamos import UnsupportedVamosFormat, UnsupportedVamosVersion, parse_native_jsonl, run_vamos_stage, select_adapter
from tr_calling_pipeline.fasta import read_fasta, sequence_sha256
from tr_calling_pipeline.schema_validation import validate

FIX=Path(__file__).parent/"fixtures"/"vamos"
def test_supported_adapter_and_unknown_rejected():
    with pytest.raises(UnsupportedVamosVersion): select_adapter("2.1.0")
    assert select_adapter("2.1.0", allow_provisional=True).capabilities.supports_bam_input
    with pytest.raises(UnsupportedVamosVersion): select_adapter("3.0.0", allow_provisional=True)
    with pytest.raises(UnsupportedVamosVersion): select_adapter(None, allow_provisional=True)
def test_lossless_known_parser():
    records=parse_native_jsonl(FIX/"supported-version"/"read.native.jsonl")
    assert [r["allele_id"] for r in records]==["allele-1","allele-2"]
    assert records[0]["repeat_count"]=="19"
    assert records[0]["motif_chain"][1]["motif"]=="CAA"
def test_malformed_rejected():
    with pytest.raises(UnsupportedVamosFormat): parse_native_jsonl(FIX/"malformed"/"malformed.jsonl")

def _terminal_config(tmp_path, *, required=False, catalog=True, opt_in=False, version="VAMOS 2.1.0"):
    fake=tmp_path/"fake vamos.py"; shutil.copy(Path(__file__).parent/"fake_vamos.py",fake); fake.chmod(0o755)
    wrapper=tmp_path/"vamos"; wrapper.write_text(f"#!/bin/sh\nFAKE_VAMOS_VERSION='{version}' exec '{fake}' \"$@\"\n"); wrapper.chmod(0o755)
    catalog_path=tmp_path/"catalog.tsv"
    if catalog: catalog_path.write_text("synthetic\n")
    locus=yaml.safe_load((Path(__file__).parents[1]/"config/loci/htt_hg38.yaml").read_text())
    locus["caller_resources"]["vamos"]["repeat_catalog"]=str(catalog_path)
    locus_path=tmp_path/"locus.yaml"; locus_path.write_text(yaml.safe_dump(locus))
    reference=tmp_path/"reference.fa"; reference.write_text(">chr4\nACGT\n"); reference_index=tmp_path/"reference.fa.fai"; reference_index.write_text("chr4\t4\t6\t4\t5\n")
    return {"run":{"case_id":"C","subject_id":"S","sample_id":"X","locus_id":"HTT"},
      "execution":{"threads":1},"inputs":{"reference_fasta":str(reference),"reference_fasta_index":str(reference_index),"assembly_records":[]},
      "locus_config":str(locus_path),"tools":{"vamos":{"executable":str(wrapper),"required":required,"allow_provisional_adapter":opt_in}}}

def _summary_schema():
    schema=json.loads((Path(__file__).parents[1]/"schemas/vamos-stage-summary.schema.json").read_text())
    # The focused validator uses local fragments, so inline the run schema here.
    run_schema=json.loads((Path(__file__).parents[1]/"schemas/vamos-run-metadata.schema.json").read_text())
    schema["properties"]["runs"]["items"]=run_schema
    return schema

def _schema(name):
    schema=json.loads((Path(__file__).parents[1]/f"schemas/{name}").read_text())
    if name=="vamos-normalized-evidence.schema.json":
        schema["properties"]["records"]["items"]=json.loads((Path(__file__).parents[1]/"schemas/normalized-caller-evidence.schema.json").read_text())
    return schema

def _validate_mode_documents(root, mode):
    directory=root/("04_vamos_read" if mode=="read" else "05_vamos_contig")
    if mode=="read":
        validate(json.loads((directory/"vamos-read.run.json").read_text()),_schema("vamos-run-metadata.schema.json"))
        validate(json.loads((directory/"vamos-read.outputs.json").read_text()),_schema("native-caller-output-registry.schema.json"))
        validate(json.loads((directory/"vamos-read.normalized.json").read_text()),_schema("vamos-normalized-evidence.schema.json"))
    else:
        validate(json.loads((directory/"stage-summary.json").read_text()),_summary_schema())
        validate(json.loads((directory/"stage-outputs.json").read_text()),_schema("native-caller-output-registry.schema.json"))
        validate(json.loads((directory/"stage-normalized.json").read_text()),_schema("vamos-normalized-evidence.schema.json"))

@pytest.mark.parametrize("kind",["unsupported","missing_catalog","missing_tool"])
def test_terminal_contig_summary_has_schema_shape(tmp_path, kind):
    config=_terminal_config(tmp_path,catalog=kind!="missing_catalog",opt_in=kind=="missing_catalog")
    if kind=="missing_tool": config["tools"]["vamos"]["executable"]=str(tmp_path/"absent")
    run_vamos_stage("04_run_vamos_contig",config,tmp_path/"run",tmp_path)
    summary=json.loads((tmp_path/"run/05_vamos_contig/stage-summary.json").read_text())
    validate(summary,_summary_schema())
    assert len(summary["runs"])==1
    expected={"unsupported":"UNSUPPORTED_VERSION","missing_catalog":"INPUT_MISSING","missing_tool":"TOOL_MISSING"}[kind]
    assert summary["runs"][0]["status"]==expected

def test_required_unverified_adapter_fails_but_preserves_metadata(tmp_path):
    config=_terminal_config(tmp_path,required=True)
    with pytest.raises(RuntimeError,match="UNSUPPORTED_VERSION"):
        run_vamos_stage("03_run_vamos_read",config,tmp_path/"run",tmp_path)
    metadata=json.loads((tmp_path/"run/04_vamos_read/vamos-read.run.json").read_text())
    assert metadata["status"]=="UNSUPPORTED_VERSION"
    assert metadata["caller_version"]=="2.1.0"
    assert "VAMOS 2.1.0" in metadata["raw_version_output"]

def test_opt_in_fake_read_and_contig_outputs_validate_and_keep_sequence_hashes(tmp_path):
    config=_terminal_config(tmp_path,opt_in=True)
    config["inputs"]["assembly_records"]=[{"record_id":"HAP1","sequence_id":"patient-hap1"},{"record_id":"HAP2","sequence_id":"patient-hap2"}]
    root=tmp_path/"run"; (root/"02_prepared_bam").mkdir(parents=True)
    (root/"02_prepared_bam/prepared.mini.bam").write_bytes(b"synthetic bam")
    (root/"02_prepared_bam/prepared.mini.bam.bai").write_bytes(b"synthetic index")
    package=root/"03_assembly_alignment"; package.mkdir()
    sequences={"patient-hap1":"ACGTACGT","patient-hap2":"ACGTCAGTCAGT"}
    (package/"patient-sequences.fasta").write_text(">patient-hap1\nACGT\nACGT\n>patient-hap2\nACGTCAGT\nCAGT\n")
    items=[]
    for index,(sequence_id,sequence) in enumerate(sequences.items(),1):
        items.append({"sequence_id":sequence_id,"source_fasta_record_id":f"HAP{index}","sequence_sha256":sequence_sha256(sequence)})
    (package/"patient-sequences.metadata.json").write_text(json.dumps({"schema_version":"1.0","sequences":items}))
    run_vamos_stage("03_run_vamos_read",config,root,tmp_path)
    run_vamos_stage("04_run_vamos_contig",config,root,tmp_path)
    _validate_mode_documents(root,"read"); _validate_mode_documents(root,"contig")
    read_records=json.loads((root/"04_vamos_read/vamos-read.normalized.json").read_text())["records"]
    assert len(read_records)==2 and all(r["assignment_state"]=="UNASSIGNED" for r in read_records)
    for sequence_id,sequence in sequences.items():
        metadata=json.loads((root/f"05_vamos_contig/{sequence_id}/run.json").read_text())
        assert metadata["sequence_sha256"]==sequence_sha256(sequence)
        assert metadata["input_fasta_sha256"] != metadata["sequence_sha256"]

def test_sequence_checksum_ignores_fasta_wrapping_but_not_bases(tmp_path):
    wrapped=tmp_path/"wrapped.fa"; unwrapped=tmp_path/"unwrapped.fa"; changed=tmp_path/"changed.fa"
    wrapped.write_text(">x\nACGT\nACGT\n"); unwrapped.write_text(">x\nACGTACGT\n"); changed.write_text(">x\nACGTACGA\n")
    digest=lambda path: sequence_sha256(read_fasta(path)[0].sequence)
    assert digest(wrapped)==digest(unwrapped)
    assert digest(wrapped)!=digest(changed)

def test_dry_run_documents_validate_for_both_modes(tmp_path):
    config=_terminal_config(tmp_path)
    root=tmp_path/"run"
    run_vamos_stage("03_run_vamos_read",config,root,tmp_path,dry_run=True)
    run_vamos_stage("04_run_vamos_contig",config,root,tmp_path,dry_run=True)
    _validate_mode_documents(root,"read"); _validate_mode_documents(root,"contig")

def test_frozen_gui_handoff_documents_validate():
    handoff=Path(__file__).parent/"fixtures/gui-handoff/task-04"
    validate(json.loads((handoff/"vamos-read.run.json").read_text()),_schema("vamos-run-metadata.schema.json"))
    validate(json.loads((handoff/"vamos-read.outputs.json").read_text()),_schema("native-caller-output-registry.schema.json"))
    validate(json.loads((handoff/"vamos-read.normalized.json").read_text()),_schema("vamos-normalized-evidence.schema.json"))
    for sequence_id in ("patient-hap1","patient-hap2"):
        directory=handoff/sequence_id
        validate(json.loads((directory/"run.json").read_text()),_schema("vamos-run-metadata.schema.json"))
        validate(json.loads((directory/"outputs.json").read_text()),_schema("native-caller-output-registry.schema.json"))
        validate(json.loads((directory/"normalized.json").read_text()),_schema("vamos-normalized-evidence.schema.json"))
    for state in ("unsupported-format","failed-call","missing-tool"):
        validate(json.loads((handoff/state/"run.json").read_text()),_schema("vamos-run-metadata.schema.json"))

def test_unsupported_native_format_preserves_and_validates_outputs(tmp_path, monkeypatch):
    config=_terminal_config(tmp_path,opt_in=True); root=tmp_path/"run"
    (root/"02_prepared_bam").mkdir(parents=True); (root/"02_prepared_bam/prepared.mini.bam").write_bytes(b"bam")
    (root/"02_prepared_bam/prepared.mini.bam.bai").write_bytes(b"index")
    monkeypatch.setenv("FAKE_VAMOS_UNSUPPORTED","1")
    run_vamos_stage("03_run_vamos_read",config,root,tmp_path)
    _validate_mode_documents(root,"read")
    metadata=json.loads((root/"04_vamos_read/vamos-read.run.json").read_text())
    normalized=json.loads((root/"04_vamos_read/vamos-read.normalized.json").read_text())
    assert metadata["status"]=="UNSUPPORTED_FORMAT"
    assert normalized["evidence_state"]=="UNSUPPORTED_FORMAT"
    assert (root/"04_vamos_read/native/vamos.jsonl").read_text()=="unsupported\n"
