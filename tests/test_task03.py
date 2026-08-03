import errno, hashlib, json, os, stat, sys
from pathlib import Path
import pytest
from tr_calling_pipeline.alignment import classify_mapping, parse_sam
from tr_calling_pipeline.bam import _link_or_copy, prepare_bam
from tr_calling_pipeline.errors import CommandExecutionError, OutputConflictError
from tr_calling_pipeline.fasta import FastaRecord, FastaValidationError, read_fasta, sequence_sha256, write_fasta
from tr_calling_pipeline.inputs import validate_reference
from tr_calling_pipeline.sequence_metadata import select_patient_sequences
from tr_calling_pipeline.tools import Tool, ToolId

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


def test_link_or_copy_conflict_overwrite_and_fallback(tmp_path, monkeypatch):
 source=tmp_path/'source'; destination=tmp_path/'destination'
 source.write_text('new'); destination.write_text('old')
 with pytest.raises(OutputConflictError,match='overwrite is disabled'):
  _link_or_copy(source,destination)
 assert destination.read_text()=='old'
 assert _link_or_copy(source,destination,overwrite=True)=='HARD_LINK'
 assert destination.read_text()=='new'
 destination.unlink()
 def cross_device(*args): raise OSError(errno.EXDEV,'cross-device link')
 monkeypatch.setattr(os,'link',cross_device)
 assert _link_or_copy(source,destination)=='COPY' and destination.read_text()=='new'


def test_link_or_copy_does_not_mask_unexpected_error(tmp_path, monkeypatch):
 source=tmp_path/'source'; destination=tmp_path/'destination'; source.write_text('data')
 def unexpected(*args): raise OSError(errno.EIO,'storage failure')
 monkeypatch.setattr(os,'link',unexpected)
 with pytest.raises(OSError) as error: _link_or_copy(source,destination)
 assert error.value.errno==errno.EIO and not destination.exists()


def fake_samtools(tmp_path):
 script=tmp_path/'fake samtools.py'
 script.write_text(f'''#!{sys.executable}
import pathlib, shutil, sys
args=sys.argv[1:]
if not args or args[0]=='--version': print('samtools 1.20'); raise SystemExit(0)
command=args[0]
if command=='quickcheck':
 raise SystemExit(0 if pathlib.Path(args[-1]).read_text().startswith('BAM:') else 2)
if command=='view': print('@HD\\tVN:1.6\\tSO:coordinate'); print('@SQ\\tSN:chrSynthetic\\tLN:32'); raise SystemExit(0)
if command=='flagstat': print('2 + 0 in total'); print('0 + 0 secondary'); print('0 + 0 supplementary'); print('2 + 0 mapped'); raise SystemExit(0)
if command=='idxstats':
 bam=pathlib.Path(args[-1]); index=pathlib.Path(str(bam)+'.bai')
 if not index.is_file() or bam.read_text().split(':',1)[1] != index.read_text().split(':',1)[1]: raise SystemExit(9)
 print('chrSynthetic\\t32\\t2\\t0'); raise SystemExit(0)
raise SystemExit(3)
''')
 script.chmod(script.stat().st_mode|stat.S_IXUSR)
 return Tool(ToolId.SAMTOOLS,'Samtools',str(script),resolved_executable=str(script),detected_version='1.20')


def test_prepare_bam_uses_exact_nonstandard_configured_index(tmp_path):
 bam=tmp_path/'bam directory'/'source.bam'; index=tmp_path/'other directory'/'chosen.index'
 bam.parent.mkdir(); index.parent.mkdir(); bam.write_text('BAM:pair-one'); index.write_text('INDEX:pair-one')
 stage=tmp_path/'output with spaces'; metadata=prepare_bam(bam,index,stage,fake_samtools(tmp_path))
 assert metadata['configured_source_index']==str(index.resolve())
 assert Path(metadata['staged_validation_index']).name=='configured-source.bam.bai'
 source_record=json.loads((stage/'execution-records/idxstats-source.json').read_text())
 assert {item['path'] for item in source_record['input_file_identities']} == {
  metadata['staged_validation_bam'], metadata['staged_validation_index']}
 prepared_record=json.loads((stage/'execution-records/idxstats-prepared.json').read_text())
 assert {item['path'] for item in prepared_record['input_file_identities']} == {
  str((stage/'prepared.mini.bam').resolve()), str((stage/'prepared.mini.bam.bai').resolve())}


def test_prepare_bam_rejects_mismatched_configured_index(tmp_path):
 bam=tmp_path/'source.bam'; index=tmp_path/'nonstandard.index'
 bam.write_text('BAM:one'); index.write_text('INDEX:different')
 with pytest.raises(CommandExecutionError): prepare_bam(bam,index,tmp_path/'stage',fake_samtools(tmp_path))


def test_prepare_bam_existing_outputs_honor_overwrite(tmp_path):
 bam=tmp_path/'source.bam'; index=tmp_path/'chosen.index'; stage=tmp_path/'stage'; stage.mkdir()
 bam.write_text('BAM:fresh'); index.write_text('INDEX:fresh')
 prepared=stage/'prepared.mini.bam'; prepared_index=stage/'prepared.mini.bam.bai'
 prepared.write_text('old bam'); prepared_index.write_text('old index')
 with pytest.raises(OutputConflictError): prepare_bam(bam,index,stage,fake_samtools(tmp_path),overwrite=False)
 assert prepared.read_text()=='old bam' and prepared_index.read_text()=='old index'
 metadata=prepare_bam(bam,index,stage,fake_samtools(tmp_path),overwrite=True)
 assert prepared.read_text()=='BAM:fresh' and prepared_index.read_text()=='INDEX:fresh'
 assert metadata['prepared_bam_sha256']==hashlib.sha256(b'BAM:fresh').hexdigest()
