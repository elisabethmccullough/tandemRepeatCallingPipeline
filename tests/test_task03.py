import hashlib, json
from pathlib import Path
import pytest
from tr_calling_pipeline.alignment import classify_mapping, parse_sam
from tr_calling_pipeline.fasta import FastaRecord, FastaValidationError, read_fasta, sequence_sha256, write_fasta
from tr_calling_pipeline.inputs import validate_reference
from tr_calling_pipeline.sequence_metadata import select_patient_sequences

FIX=Path(__file__).parent/'fixtures/inputs'

def configured(*ids):
 return [{'record_id':rid,'sequence_id':f'patient-{rid.lower()}','display_label':rid,'sequence_role':'PATIENT_CONSENSUS' if rid=='CONSENSUS' else 'PATIENT_HAPLOTYPE'} for rid in ids]

def test_fasta_validation_selection_hash_and_output(tmp_path):
 records=read_fasta(FIX/'assembly-two-haplotypes.fasta'); assert [r.identifier for r in records]==['HAP1','HAP2','UNRELATED']
 selected,metadata=select_patient_sequences(FIX/'assembly-two-haplotypes.fasta',configured('HAP2','HAP1'))
 assert [r.identifier for r in selected]==['patient-hap2','patient-hap1'] and metadata[0].source_fasta_record_id=='HAP2'
 out=tmp_path/'patient fasta.fa'; write_fasta(out,selected,width=5)
 assert out.read_text().startswith('>patient-hap2\nTTGCA\n') and read_fasta(out)[0].sequence=='TTGCAACGTTAA'
 assert sequence_sha256('acgt')==hashlib.sha256(b'ACGT').hexdigest()

def test_consensus_lowercase_iupac():
 selected,metadata=select_patient_sequences(FIX/'assembly-consensus-only.fasta',configured('CONSENSUS'))
 assert selected[0].sequence=='acgtryN' and metadata[0].mapping_status=='NOT_ALIGNED'

@pytest.mark.parametrize('name', ['assembly-duplicate-ids.fasta','assembly-empty-record.fasta','assembly-invalid-character.fasta'])
def test_bad_fasta(name):
 with pytest.raises(FastaValidationError): read_fasta(FIX/name)

def test_missing_and_duplicate_selection():
 with pytest.raises(FastaValidationError,match='not found'): select_patient_sequences(FIX/'assembly-two-haplotypes.fasta',configured('NOPE'))
 bad=configured('HAP1','HAP2'); bad[1]['sequence_id']=bad[0]['sequence_id']
 with pytest.raises(FastaValidationError,match='sequence_id'): select_patient_sequences(FIX/'assembly-two-haplotypes.fasta',bad)

def test_reference_validation():
 report=validate_reference(FIX/'reference-local.fasta',FIX/'reference-local.fasta.fai',scope='LOCAL_LOCUS_REFERENCE',required_contig='chrSynthetic')
 assert report['reference_scope']=='LOCAL_LOCUS_REFERENCE' and len(report['reference_fasta_sha256'])==64
 with pytest.raises(ValueError,match='absent'): validate_reference(FIX/'reference-local.fasta',FIX/'reference-local.fasta.fai',required_contig='chrMissing')

def test_alignment_orientation_coordinates_and_states():
 lines=['q\t0\tchr1\t5\t42\t4M2D3M\t*\t0\t0\tACGTAAA\t*','r\t16\tchr1\t20\t60\t7M\t*\t0\t0\tAAAAAAA\t*','u\t4\t*\t0\t0\t*\t*\t0\t0\tNNNN\t*','s\t256\tchr1\t3\t10\t4M\t*\t0\t0\tAAAA\t*']
 alns=parse_sam(lines)
 q=classify_mapping('q','Q','0'*64,7,alns); assert (q.alignment_start,q.alignment_end,q.mapped_base_count)==(4,13,7) and q.reverse_complement_required is False
 r=classify_mapping('r','R','1'*64,7,alns); assert r.reference_alignment_strand=='REVERSE' and r.reverse_complement_required is True
 assert classify_mapping('u','U','2'*64,4,alns).mapping_status=='UNMAPPED'
 assert classify_mapping('s','S','3'*64,4,alns).mapping_status=='SECONDARY_ONLY'

def test_new_json_fixtures_parse():
 for path in (Path(__file__).parent/'fixtures/gui-handoff/task-03').glob('*.json'): json.loads(path.read_text())
