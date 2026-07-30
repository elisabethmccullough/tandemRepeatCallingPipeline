import csv, json
from tr_calling_pipeline.manifest import create_manifest, write_manifest
from tr_calling_pipeline.normalize import NORMALIZED_COLUMNS, write_normalized_outputs

def test_empty_normalized_outputs(tmp_path):
    paths=write_normalized_outputs(tmp_path)
    with paths['tsv'].open() as handle: assert next(csv.reader(handle,delimiter='\t'))==NORMALIZED_COLUMNS
    assert json.loads(paths['json'].read_text())==[]
    assert paths['warnings'].exists()

def test_manifest_creation(tmp_path):
    source=tmp_path/'a.fa'; source.write_text('>a\nA\n')
    config={'run':{'sample_id':'S1','locus_id':'L1'},'inputs':{'assembly_fasta':str(source)}}
    manifest=create_manifest(config,tmp_path/'config.yaml')
    path=write_manifest(manifest,tmp_path/'run_manifest.json')
    loaded=json.loads(path.read_text())
    assert loaded['sample_id']=='S1'; assert loaded['input_checksums']['assembly_fasta']
