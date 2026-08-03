from pathlib import Path
import pytest
from tr_calling_pipeline.last_alignment import classify_version as lc, select_adapter as la, validate_maf, write_single_record_fasta, LastVersionClassification, UnsupportedLastVersion
from tr_calling_pipeline.callers.tandem_genotypes import classify_version as tc, select_adapter as ta, parse_native_tsv, normalize_record, TandemGenotypesVersionClassification, UnsupportedTandemGenotypesVersion
from tr_calling_pipeline.caller_outputs import NativeCallerOutput

FIX=Path(__file__).parent/'fixtures/tandem-genotypes'
def test_adapters_are_explicitly_gated():
 assert lc('1450') is LastVersionClassification.PROVISIONAL_DEVELOPMENT
 with pytest.raises(UnsupportedLastVersion): la('1450')
 assert la('1450',allow_provisional=True).lastdb_plan('lastdb',Path('db'),Path('in.fa')).argv == ('lastdb','db','in.fa')
 assert tc('0.1.0') is TandemGenotypesVersionClassification.PROVISIONAL_DEVELOPMENT
 with pytest.raises(UnsupportedTandemGenotypesVersion): ta('0.1.0')
 assert ta('0.1.0',allow_provisional=True).plan('tg',alignment=Path('a.maf'),repeat_definition=Path('r'),output=Path('o')).argv == ('tg','a.maf','r','o')
def test_maf_and_deterministic_fasta(tmp_path):
 assert 'patient-hap1' in validate_maf(FIX/'last/valid-alignment.maf','patient-hap1')
 p=tmp_path/'space here'/'patient-hap1.fasta'; digest=write_single_record_fasta(p,'patient-hap1','CAG'*30)
 assert p.read_text().startswith('>patient-hap1\n') and len(digest)==64
@pytest.mark.parametrize('name',['malformed-alignment.maf','unsupported-alignment.txt'])
def test_invalid_maf(name):
 with pytest.raises(Exception): validate_maf(FIX/'last'/name,'patient-hap1')
def test_lossless_direct_normalization():
 p=FIX/'native/patient-hap1.native.tsv'; records=parse_native_tsv(p)
 assert records[0].raw_fields['repeat_count']=='020' and records[0].raw_fields['caller_extra']=='00123'
 source=NativeCallerOutput.from_path(p,file_id='native-h1',caller='TANDEM_GENOTYPES',caller_version='0.1.0',analysis_source='ASSEMBLED_CONTIG',producer_command_id='cmd')
 c={'run':{'case_id':'case','subject_id':'subject','sample_id':'sample','locus_id':'HTT'}}
 n=normalize_record(records[0],config=c,caller_version='0.1.0',source=source,associated_sequence_id='patient-hap1')
 assert n['caller']=='TANDEM_GENOTYPES' and n['assignment_state']=='DIRECT_SEQUENCE_ASSOCIATION'
 assert n['associated_sequence_id']=='patient-hap1' and n['start'] is None and n['raw_fields']['repeat_count']=='020'
 with pytest.raises(ValueError): normalize_record(records[0],config=c,caller_version='0.1.0',source=source,associated_sequence_id='patient-hap2')
