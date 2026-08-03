from pathlib import Path
from copy import deepcopy
import pytest
import yaml
from tr_calling_pipeline.config import ConfigurationError, load_config, validate_identifier
from tr_calling_pipeline.manifest import create_manifest

BASE = yaml.safe_load((Path(__file__).parents[1] / "config/example.yaml").read_text())

def test_load_valid_configuration(tmp_path):
    path=tmp_path/'config.yaml'; path.write_text(yaml.safe_dump(BASE)); config=load_config(path)
    assert config['run']['sample_id']=='HG00438'; assert Path(config['inputs']['assembly_fasta']).is_absolute()

def test_paths_resolve_from_configuration_directory(tmp_path, monkeypatch):
    config_dir = tmp_path / "configuration"; config_dir.mkdir()
    document = deepcopy(BASE)
    for key in ('assembly_fasta','mini_bam','mini_bam_index','reference_fasta','reference_fasta_index'):
        document['inputs'][key] = f'inputs/{key}'
        target=config_dir/document['inputs'][key]; target.parent.mkdir(exist_ok=True); target.write_text('fixture')
    document['locus_config']='loci/locus.yaml'; locus=config_dir/'loci/locus.yaml'; locus.parent.mkdir(); locus.write_text('fixture')
    document['run']['output_root']='results'; document['case_package']['package_root']='packages'
    path=config_dir/'run.yaml'; path.write_text(yaml.safe_dump(document))
    elsewhere=tmp_path/'elsewhere'; elsewhere.mkdir(); monkeypatch.chdir(elsewhere)
    config=load_config(path,check_inputs=True)
    assert config['inputs']['assembly_fasta']==str((config_dir/'inputs/assembly_fasta').resolve())
    assert config['locus_config']==str(locus.resolve())
    manifest=create_manifest(config,path)
    assert manifest['resolved_paths']['output_root']==str((config_dir/'results').resolve())

def test_reject_missing_required_field(tmp_path):
    broken=deepcopy(BASE); del broken['run']['locus_id']; path=tmp_path/'bad.yaml'; path.write_text(yaml.safe_dump(broken))
    with pytest.raises(ConfigurationError, match='locus_id'): load_config(path)

@pytest.mark.parametrize('value',['sample_1','HTT-1','A.b'])
def test_safe_identifiers(value): assert validate_identifier(value)==value
@pytest.mark.parametrize('value',['../bad','has space',''])
def test_unsafe_identifiers(value):
    with pytest.raises(ConfigurationError): validate_identifier(value)
