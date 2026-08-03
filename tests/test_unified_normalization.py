import hashlib, json
from pathlib import Path
from tr_calling_pipeline.normalization import normalize_caller_artifacts
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
