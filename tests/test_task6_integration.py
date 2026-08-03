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
