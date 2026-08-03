import hashlib,json,shutil
from pathlib import Path
import yaml
from tr_calling_pipeline.last_alignment import prepare_tandem_genotypes_stage
from tr_calling_pipeline.callers.tandem_genotypes import run_tandem_genotypes_stage

def test_two_sequence_fake_execution_and_stale_isolation(tmp_path):
 bindir=tmp_path/'tools with spaces'; bindir.mkdir()
 for name in ('lastdb','lastal','tandem-genotypes'):
  p=bindir/name; shutil.copy(Path(__file__).parent/'fake_task6_tool.py',p); p.chmod(0o755)
 root=tmp_path/'run'; source=root/'03_assembly_alignment'; source.mkdir(parents=True)
 seqs={'patient-hap1':'CAGCAG','patient-hap2':'CAACAA'}
 (source/'patient-sequences.fasta').write_text(''.join(f'>{k}\n{v}\n' for k,v in seqs.items()))
 meta={'schema_version':'1.0','sequences':[{'sequence_id':k,'source_fasta_record_id':'source-'+k,'sequence_sha256':hashlib.sha256(v.encode()).hexdigest()} for k,v in seqs.items()]}
 (source/'patient-sequences.metadata.json').write_text(json.dumps(meta))
 ref=tmp_path/'ref.fa'; ref.write_text('>chr4\nCAG\n'); (tmp_path/'ref.fa.fai').write_text('chr4\t3\t6\t3\t4\n'); repeat=tmp_path/'repeat.txt'; repeat.write_text('synthetic\n')
 locus=yaml.safe_load(Path('config/loci/htt_hg38.yaml').read_text()); locus['caller_resources']['tandem_genotypes'].update(repeat_definition=str(repeat),allow_provisional_last_adapter=True,allow_provisional_tandem_genotypes_adapter=True)
 lp=tmp_path/'locus.yaml'; lp.write_text(yaml.safe_dump(locus))
 config={'locus_config':str(lp),'inputs':{'reference_fasta':str(ref),'reference_fasta_index':str(tmp_path/'ref.fa.fai')},'execution':{'threads':1},'run':{'case_id':'case','subject_id':'subject','sample_id':'sample','locus_id':'HTT'},'tools':{'lastdb':{'executable':str(bindir/'lastdb'),'required':True},'lastal':{'executable':str(bindir/'lastal'),'required':True},'tandem_genotypes':{'executable':str(bindir/'tandem-genotypes'),'required':True}}}
 prepare_tandem_genotypes_stage(config,root,tmp_path)
 stale=root/'08_tandem_genotypes'/'patient-hap1'/'native'/'stale.txt'
 run_tandem_genotypes_stage(config,root,tmp_path)
 assert not stale.exists()
 normalized=json.loads((root/'08_tandem_genotypes/stage-normalized.json').read_text())['records']
 assert {r['associated_sequence_id'] for r in normalized}==set(seqs)
 assert all(r['assignment_state']=='DIRECT_SEQUENCE_ASSOCIATION' for r in normalized)
 assert len(json.loads((root/'07_tandem_genotypes_preparation/alignment-inputs.json').read_text())['alignments'])==2
 assert json.loads((root/'08_tandem_genotypes/stage-summary.json').read_text())['status']=='SUCCEEDED'

def _dry_config(tmp_path, *, repeat=True, lastdb=True, lastal=True, provisional=True):
 bindir=tmp_path/'dry tools'; bindir.mkdir()
 paths={}
 for name,available in (('lastdb',lastdb),('lastal',lastal)):
  p=bindir/name
  if available: shutil.copy(Path(__file__).parent/'fake_task6_tool.py',p); p.chmod(0o755)
  paths[name]=p
 root=tmp_path/'dry-run'; source=root/'03_assembly_alignment'; source.mkdir(parents=True)
 seq='CAG'; (source/'patient-sequences.fasta').write_text('>patient-hap1\n'+seq+'\n')
 (source/'patient-sequences.metadata.json').write_text(json.dumps({'sequences':[{'sequence_id':'patient-hap1','source_fasta_record_id':'source-1','sequence_sha256':hashlib.sha256(seq.encode()).hexdigest()}]}))
 ref=tmp_path/'dry-ref.fa'; ref.write_text('>chr4\nCAG\n'); definition=tmp_path/'dry-repeat.txt'
 if repeat: definition.write_text('synthetic\n')
 locus=yaml.safe_load(Path('config/loci/htt_hg38.yaml').read_text()); locus['caller_resources']['tandem_genotypes'].update(repeat_definition=str(definition) if repeat else None,allow_provisional_last_adapter=provisional)
 lp=tmp_path/'dry-locus.yaml'; lp.write_text(yaml.safe_dump(locus))
 config={'locus_config':str(lp),'inputs':{'reference_fasta':str(ref),'reference_fasta_index':str(tmp_path/'none.fai')},'tools':{'lastdb':{'executable':str(paths['lastdb']),'required':True},'lastal':{'executable':str(paths['lastal']),'required':True}}}
 return config,root

import pytest
@pytest.mark.parametrize('options',[{'repeat':False},{'lastdb':False},{'lastal':False},{'provisional':False},{}])
def test_prepare_dry_run_always_remains_dry_run(tmp_path, options):
 config,root=_dry_config(tmp_path,**options)
 records=prepare_tandem_genotypes_stage(config,root,tmp_path,dry_run=True)
 summary=json.loads((root/'07_tandem_genotypes_preparation/preparation-summary.json').read_text())
 assert summary['status']=='DRY_RUN' and all(r['status']=='DRY_RUN' for r in records)
 assert not list((root/'07_tandem_genotypes_preparation').glob('*/input/*.fasta'))
 assert summary['warnings']

def test_prepare_dry_run_unsupported_version_remains_dry_run(tmp_path,monkeypatch):
 config,root=_dry_config(tmp_path); monkeypatch.setenv('FAKE_VERSION','LAST unknown')
 records=prepare_tandem_genotypes_stage(config,root,tmp_path,dry_run=True)
 assert records[0]['status']=='DRY_RUN' and 'version' in records[0]['warnings'][0].lower()

@pytest.mark.parametrize('preparations,expected',[([], 'NOT_COMPUTED'),([{'status':'FAILED'}],'NOT_COMPUTED'),([{'status':'INPUT_MISSING'}],'INPUT_MISSING'),([{'status':'TOOL_MISSING'}],'TOOL_MISSING'),([{'status':'UNSUPPORTED_VERSION'}],'UNSUPPORTED_VERSION'),([{'status':'UNSUPPORTED_FORMAT'}],'UNSUPPORTED_FORMAT')])
def test_run_with_no_eligible_preparation_is_terminal(tmp_path,preparations,expected):
 config,root=_dry_config(tmp_path); prep=root/'07_tandem_genotypes_preparation'; prep.mkdir(parents=True); (prep/'alignment-metadata.json').write_text(json.dumps({'records':preparations}))
 config['tools']['tandem_genotypes']={'executable':str(tmp_path/'missing-tg'),'required':False}
 assert run_tandem_genotypes_stage(config,root,tmp_path)==[]
 summary=json.loads((root/'08_tandem_genotypes/stage-summary.json').read_text()); normalized=json.loads((root/'08_tandem_genotypes/stage-normalized.json').read_text())
 assert summary['status']==expected and summary['runs']==[] and normalized['records']==[]
 assert expected!='SUCCEEDED'
