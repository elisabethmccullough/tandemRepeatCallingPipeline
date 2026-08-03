import json, shutil
from dataclasses import asdict
from pathlib import Path
import pytest, yaml
from tr_calling_pipeline.callers.straglr import (UnsupportedStraglrFormat, UnsupportedStraglrVersion,
    classify_version, parse_native_tsv, run_straglr_stage, select_adapter, StraglrVersionClassification)
from tr_calling_pipeline.schema_validation import validate
from tr_calling_pipeline.provenance import file_identity
from tr_calling_pipeline.stages import ResumeReason, resume_eligibility

FIX=Path(__file__).parent/'fixtures/straglr'
def test_adapter_gating_and_classification():
    assert classify_version('1.5.0') is StraglrVersionClassification.PROVISIONAL_DEVELOPMENT
    assert classify_version(None) is StraglrVersionClassification.UNDETERMINED
    assert classify_version('9.0.0') is StraglrVersionClassification.UNSUPPORTED
    with pytest.raises(UnsupportedStraglrVersion): select_adapter('1.5.0')
    assert select_adapter('1.5.0',allow_provisional=True).capabilities.native_coordinate_convention=='zero_based_half_open'
    with pytest.raises(UnsupportedStraglrVersion): select_adapter('9.0.0',allow_provisional=True)
def test_lossless_parser_and_malformed():
    rows=parse_native_tsv(FIX/'supported-version/two-alleles.native.tsv')
    assert [r.native_allele_identifier for r in rows]==['allele-1','allele-2']
    assert rows[0].raw_fields['repeat_count']=='004' and rows[0].raw_fields['caller_specific']=='0007'
    with pytest.raises(UnsupportedStraglrFormat): parse_native_tsv(FIX/'malformed/malformed.native.tsv')

def config(tmp_path, *, required=False, opt_in=True, catalog=True):
    tmp_path.mkdir(parents=True,exist_ok=True)
    fake=tmp_path/'fake straglr.py'; shutil.copy(Path(__file__).parent/'fake_straglr.py',fake); fake.chmod(0o755)
    cat=tmp_path/'repeat catalog.tsv'
    if catalog: cat.write_text('synthetic catalog\n')
    locus=yaml.safe_load((Path(__file__).parents[1]/'config/loci/htt_hg38.yaml').read_text())
    locus['caller_resources']['straglr']['repeat_catalog']=cat.name
    locus['caller_resources']['straglr']['allow_provisional_adapter']=opt_in
    locus['caller_resources']['straglr']['additional_arguments']=['--synthetic-flag']
    lp=tmp_path/'locus config.yaml'; lp.write_text(yaml.safe_dump(locus))
    ref=tmp_path/'reference genome.fa'; ref.write_text('>chrSynthetic\nACGT\n'); (tmp_path/'reference genome.fa.fai').write_text('chrSynthetic\t4\t14\t4\t5\n')
    return {'run':{'case_id':'CASE','subject_id':'SUBJECT','sample_id':'SAMPLE','locus_id':'SYNTH_LOCUS'},'execution':{'threads':2},
      'inputs':{'reference_fasta':str(ref),'reference_fasta_index':str(ref)+'.fai'},'locus_config':str(lp),
      'tools':{'straglr':{'executable':str(fake),'required':required}}}
def prepared(root):
    p=root/'02_prepared_bam'; p.mkdir(parents=True); (p/'prepared.mini.bam').write_bytes(b'bam'); (p/'prepared.mini.bam.bai').write_bytes(b'index')
def test_execution_normalization_registration_and_immutability(tmp_path,monkeypatch):
    cfg=config(tmp_path); root=tmp_path/'output with spaces'; prepared(root)
    catalog=tmp_path/'repeat catalog.tsv'; before=catalog.read_bytes()
    monkeypatch.setenv('FAKE_STRAGLR_MULTIPLE','1')
    run_straglr_stage(cfg,root,tmp_path)
    doc=json.loads((root/'06_straglr/straglr.normalized.json').read_text()); outputs=json.loads((root/'06_straglr/straglr.outputs.json').read_text())
    assert len(doc['records'])==2 and len(outputs['outputs'])==2
    assert all(x['caller']=='STRAGLR' and x['analysis_source']=='RAW_READS' for x in outputs['outputs'])
    assert all(x['assignment_state']=='UNASSIGNED' and x['associated_sequence_id'] is None for x in doc['records'])
    assert doc['records'][0]['raw_fields']['repeat_count']=='004'
    assert catalog.read_bytes()==before
    native=root/'06_straglr/native/straglr.tsv'; digest=outputs['outputs'][1 if outputs['outputs'][1]['path'].endswith('straglr.tsv') else 0]['sha256']
    import hashlib
    assert hashlib.sha256(native.read_bytes()).hexdigest()==digest
    assert str(root/'02_prepared_bam/prepared.mini.bam') in json.loads((root/'06_straglr/execution-records/straglr-read.json').read_text())['argv']
    assert '--synthetic-flag' in json.loads((root/'06_straglr/execution-records/straglr-read.json').read_text())['argv']
def test_terminal_states_and_dry_run(tmp_path):
    root=tmp_path/'root'; prepared(root); cfg=config(tmp_path,opt_in=False)
    run_straglr_stage(cfg,root,tmp_path); assert json.loads((root/'06_straglr/straglr.run.json').read_text())['status']=='UNSUPPORTED_VERSION'
    root2=tmp_path/'root2'; prepared(root2); cfg2=config(tmp_path/'second',catalog=False); run_straglr_stage(cfg2,root2,tmp_path/'second')
    assert json.loads((root2/'06_straglr/straglr.run.json').read_text())['status']=='INPUT_MISSING'
    root3=tmp_path/'root3'; run_straglr_stage(cfg,root3,tmp_path,dry_run=True)
    assert json.loads((root3/'06_straglr/straglr.normalized.json').read_text())['evidence_state']=='NOT_COMPUTED'
def test_unsupported_format_keeps_native(tmp_path,monkeypatch):
    cfg=config(tmp_path); root=tmp_path/'root'; prepared(root); monkeypatch.setenv('FAKE_STRAGLR_UNSUPPORTED','1')
    run_straglr_stage(cfg,root,tmp_path)
    assert (root/'06_straglr/native/straglr.tsv').read_text()=='unknown\tcolumns\nx\ty\n'
    assert json.loads((root/'06_straglr/straglr.normalized.json').read_text())['evidence_state']=='UNSUPPORTED_FORMAT'

def test_native_output_set_isolated_on_overwrite(tmp_path,monkeypatch):
    cfg=config(tmp_path); root=tmp_path/'root'; prepared(root)
    run_straglr_stage(cfg,root,tmp_path)
    stale=root/'06_straglr/native/stale.txt'; stale.write_text('old run')
    with pytest.raises(FileExistsError): run_straglr_stage(cfg,root,tmp_path)
    monkeypatch.setenv('FAKE_STRAGLR_MULTIPLE','1')
    run_straglr_stage(cfg,root,tmp_path,overwrite=True)
    registry=json.loads((root/'06_straglr/straglr.outputs.json').read_text())['outputs']
    assert {Path(item['path']).name for item in registry}=={'straglr.tsv','straglr.details.txt'}
    assert not stale.exists()

def test_failed_overwrite_keeps_old_native_set(tmp_path,monkeypatch):
    cfg=config(tmp_path); root=tmp_path/'root'; prepared(root); run_straglr_stage(cfg,root,tmp_path)
    old=(root/'06_straglr/native/straglr.tsv').read_bytes()
    monkeypatch.setenv('FAKE_STRAGLR_FAIL','1')
    with pytest.raises(Exception): run_straglr_stage(cfg,root,tmp_path,overwrite=True)
    assert (root/'06_straglr/native/straglr.tsv').read_bytes()==old
    assert not list((root/'06_straglr').glob('.native-work-*'))

def test_run_schema_rejects_ignored_tool_runtime_settings():
    schema=json.loads((Path(__file__).parents[1]/'schemas/run-config.schema.json').read_text())
    straglr=schema['properties']['tools']['properties']['straglr']
    assert 'allow_provisional_adapter' not in straglr['properties']
    assert 'additional_arguments' not in straglr['properties']
    assert straglr['additionalProperties'] is False

@pytest.mark.parametrize('field,value',[('allow_provisional_adapter',False),('additional_arguments',['--changed'])])
def test_locus_runtime_setting_change_invalidates_resume(tmp_path,field,value):
    cfg=config(tmp_path); locus=Path(cfg['locus_config']); before=asdict(file_identity(locus))
    prior={'record_schema_version':'1.0','status':'SUCCEEDED','configuration_digest':'digest',
      'input_file_identities':[before],'output_file_identities':[],'tool_identities':[]}
    document=yaml.safe_load(locus.read_text()); document['caller_resources']['straglr'][field]=value
    locus.write_text(yaml.safe_dump(document))
    reason=resume_eligibility(prior,'digest',[asdict(file_identity(locus))],[],[])
    assert reason is ResumeReason.INPUT_CHANGED
