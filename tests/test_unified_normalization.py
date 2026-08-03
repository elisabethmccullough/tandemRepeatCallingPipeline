import hashlib, json
from pathlib import Path
import pytest
from tr_calling_pipeline.normalization import NormalizationError, normalize_caller_artifacts, run_normalization_stage
from tr_calling_pipeline.schema_validation import validate


def test_all_callers_unavailable_is_valid_and_deterministic(tmp_path):
    seq='ACGT'; fasta=tmp_path/'patient.fa'; fasta.write_text('>patient-hap1\n'+seq+'\n')
    metadata={'schema_version':'1.0','sequences':[{'record_schema_version':'1.0','sequence_id':'patient-hap1','sequence_role':'PATIENT_HAPLOTYPE','display_label':'Patient sequence 1','source_fasta_record_id':'patient-hap1','sequence_sha256':hashlib.sha256(seq.encode()).hexdigest(),'sequence_length':4,'original_orientation':'AS_STORED','reference_alignment_strand':'NOT_ALIGNED','display_orientation':'UNRESOLVED','reverse_complement_required':None,'source_coordinates':None,'mapping_status':'NOT_ALIGNED','warnings':[]}]}
    meta=tmp_path/'metadata.json'; meta.write_text(json.dumps(metadata))
    args=dict(case_id='C',subject_id='P',sample_id='S',locus_id='L',patient_metadata_path=meta,patient_fasta_path=fasta,artifact_groups=[],created_utc='2026-01-01T00:00:00Z')
    first=normalize_caller_artifacts(**args); second=normalize_caller_artifacts(**args)
    assert first == second
    assert first[0]['records'] == []
    assert [x['evidence_state'] for x in first[0]['caller_summaries']] == ['NOT_COMPUTED']*4
    root=Path(__file__).parents[1]
    for document,name in zip(first,('unified-normalized-evidence.schema.json','unified-evidence-summary.schema.json','unified-source-registry.schema.json','unified-validation-report.schema.json')):
        validate(document,json.loads((root/'schemas'/name).read_text()))


def test_gui_fixture_has_stable_caller_distinctions():
    root=Path(__file__).parent/'fixtures/gui-handoff/task-07'
    package=json.loads((root/'normalized-evidence.json').read_text())
    triples={(r['caller'],r['analysis_source'],r['assignment_state'],r['associated_sequence_id']) for r in package['records']}
    assert ('VAMOS','RAW_READS','UNASSIGNED',None) in triples
    assert ('STRAGLR','RAW_READS','UNASSIGNED',None) in triples
    for sequence in ('patient-hap1','patient-hap2'):
        assert ('VAMOS','ASSEMBLED_CONTIG','DIRECT_SEQUENCE_ASSOCIATION',sequence) in triples
        assert ('TANDEM_GENOTYPES','ASSEMBLED_CONTIG','DIRECT_SEQUENCE_ASSOCIATION',sequence) in triples


def _stage_inputs(tmp_path):
    root=tmp_path/'run'; patient=root/'03_assembly_alignment'; patient.mkdir(parents=True)
    sequence='ACGT'; (patient/'patient-sequences.fasta').write_text(f'>patient-hap1\n{sequence}\n')
    metadata={'schema_version':'1.0','sequences':[{'record_schema_version':'1.0','sequence_id':'patient-hap1','sequence_role':'PATIENT_HAPLOTYPE','display_label':'Patient sequence 1','source_fasta_record_id':'patient-hap1','sequence_sha256':hashlib.sha256(sequence.encode()).hexdigest(),'sequence_length':4,'original_orientation':'AS_STORED','reference_alignment_strand':'NOT_ALIGNED','display_orientation':'UNRESOLVED','reverse_complement_required':None,'source_coordinates':None,'mapping_status':'NOT_ALIGNED','warnings':[]}]}
    (patient/'patient-sequences.metadata.json').write_text(json.dumps(metadata))
    locus=tmp_path/'locus.yaml'; locus.write_text('schema_version: "1.0"\n')
    config={'run':{'case_id':'C','subject_id':'P','sample_id':'S','locus_id':'L'},'locus_config':str(locus)}
    return root,config


def test_stage_all_optional_groups_absent_become_not_computed(tmp_path):
    root,config=_stage_inputs(tmp_path)
    paths=run_normalization_stage(config,root)
    package=json.loads(paths['normalized-evidence.json'].read_text())
    assert [item['evidence_state'] for item in package['caller_summaries']] == ['NOT_COMPUTED']*4


def test_stage_dry_run_reports_absent_and_incomplete_without_package(tmp_path):
    root,config=_stage_inputs(tmp_path)
    partial=root/'06_straglr'; partial.mkdir(parents=True); (partial/'straglr.outputs.json').write_text('{}')
    paths=run_normalization_stage(config,root,dry_run=True)
    plan=json.loads(paths['dry-run-validation-report.json'].read_text())
    statuses={(item['caller'],item['analysis_source']):item['status'] for item in plan['caller_artifact_groups']}
    assert statuses[('STRAGLR','RAW_READS')] == 'INCOMPLETE'
    assert statuses[('TANDEM_GENOTYPES','ASSEMBLED_CONTIG')] == 'ABSENT'
    assert not plan['normalization_performed'] and not plan['valid_for_execution']
    assert not (root/'09_normalized_evidence/normalized-evidence.json').exists()
    assert not (root/'09_normalized_evidence/source-registry.json').exists()


@pytest.mark.parametrize('present_directories,expected_absent', [
    (('04_vamos_read','06_straglr'), {'VAMOS:ASSEMBLED_CONTIG','TANDEM_GENOTYPES:ASSEMBLED_CONTIG'}),
    (('05_vamos_contig','08_tandem_genotypes'), {'VAMOS:RAW_READS','STRAGLR:RAW_READS'}),
    (('04_vamos_read','05_vamos_contig','08_tandem_genotypes'), {'STRAGLR:RAW_READS'}),
    (('04_vamos_read','05_vamos_contig','06_straglr'), {'TANDEM_GENOTYPES:ASSEMBLED_CONTIG'}),
])
def test_dry_run_optional_caller_combinations(tmp_path,present_directories,expected_absent):
    root,config=_stage_inputs(tmp_path)
    names={'04_vamos_read':('vamos-read.normalized.json','vamos-read.outputs.json','vamos-read.run.json'),
           '05_vamos_contig':('stage-normalized.json','stage-outputs.json','stage-summary.json'),
           '06_straglr':('straglr.normalized.json','straglr.outputs.json','straglr.run.json'),
           '08_tandem_genotypes':('stage-normalized.json','stage-outputs.json','stage-summary.json')}
    for directory in present_directories:
        path=root/directory; path.mkdir(parents=True)
        for name in names[directory]: (path/name).write_text('{}')
    plan=json.loads(run_normalization_stage(config,root,dry_run=True)['dry-run-validation-report.json'].read_text())
    absent={f'{item["caller"]}:{item["analysis_source"]}' for item in plan['caller_artifact_groups'] if item['status']=='ABSENT'}
    assert absent == expected_absent
    assert all(item['status'] != 'INCOMPLETE' for item in plan['caller_artifact_groups'])


@pytest.mark.parametrize('required', ['patient-sequences.fasta','patient-sequences.metadata.json'])
def test_stage_checks_required_patient_inputs(tmp_path,required):
    root,config=_stage_inputs(tmp_path); (root/'03_assembly_alignment'/required).unlink()
    with pytest.raises(FileNotFoundError,match='required Stage 08 input'):
        run_normalization_stage(config,root)
    plan=json.loads(run_normalization_stage(config,root,dry_run=True)['dry-run-validation-report.json'].read_text())
    assert not plan['valid_for_execution']
    assert any(required in item['path'] and not item['exists'] for item in plan['required_inputs'])


def test_present_corrupt_and_partial_groups_fail(tmp_path):
    root,config=_stage_inputs(tmp_path)
    group=root/'04_vamos_read'; group.mkdir(parents=True)
    for name in ('vamos-read.normalized.json','vamos-read.outputs.json','vamos-read.run.json'):
        (group/name).write_text('{not-json')
    with pytest.raises(NormalizationError): run_normalization_stage(config,root)
    assert not (root/'09_normalized_evidence').exists()
    reports=list((root/'09_normalized_evidence_failures').glob('*/validation-report.json'))
    assert len(reports)==1


def test_failed_overwrite_preserves_completed_package(tmp_path):
    root,config=_stage_inputs(tmp_path); run_normalization_stage(config,root)
    completed=(root/'09_normalized_evidence/normalized-evidence.json').read_bytes()
    group=root/'04_vamos_read'; group.mkdir(parents=True)
    for name in ('vamos-read.normalized.json','vamos-read.outputs.json','vamos-read.run.json'):
        (group/name).write_text('{not-json')
    with pytest.raises(NormalizationError): run_normalization_stage(config,root,overwrite=True)
    assert (root/'09_normalized_evidence/normalized-evidence.json').read_bytes()==completed
