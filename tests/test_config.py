from pathlib import Path
import pytest
import yaml
from tr_calling_pipeline.config import ConfigurationError, load_config, validate_identifier

BASE = {"run":{"sample_id":"S1","locus_id":"L1","output_root":"outputs"},"inputs":{"assembly_fasta":"a.fa","mini_bam":"a.bam","mini_bam_index":"a.bam.bai","reference_fasta":"r.fa","reference_fasta_index":"r.fa.fai"},"locus_config":"locus.yaml","execution":{"threads":1,"overwrite":False}}

def test_load_valid_configuration(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path); path=tmp_path/'config.yaml'; path.write_text(yaml.safe_dump(BASE))
    config=load_config(path)
    assert config['run']['sample_id']=='S1'; assert Path(config['inputs']['assembly_fasta']).is_absolute()

def test_reject_missing_required_field(tmp_path):
    broken={**BASE,"run":{"sample_id":"S1","output_root":"outputs"}}; path=tmp_path/'bad.yaml'; path.write_text(yaml.safe_dump(broken))
    with pytest.raises(ConfigurationError, match='run.locus_id'): load_config(path)

@pytest.mark.parametrize('value',['sample_1','HTT-1','A.b'])
def test_safe_identifiers(value): assert validate_identifier(value)==value
@pytest.mark.parametrize('value',['../bad','has space',''])
def test_unsafe_identifiers(value):
    with pytest.raises(ConfigurationError): validate_identifier(value)
