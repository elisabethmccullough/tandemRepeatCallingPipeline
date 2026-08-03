from pathlib import Path
from copy import deepcopy
import pytest
import yaml
from tr_calling_pipeline.config import ConfigurationError, load_config, validate_identifier
from tr_calling_pipeline.manifest import create_manifest

BASE = {"run":{"sample_id":"S1","locus_id":"L1","output_root":"outputs"},"inputs":{"assembly_fasta":"a.fa","mini_bam":"a.bam","mini_bam_index":"a.bam.bai","reference_fasta":"r.fa","reference_fasta_index":"r.fa.fai"},"locus_config":"locus.yaml","execution":{"threads":1,"overwrite":False}}

def test_load_valid_configuration(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path); path=tmp_path/'config.yaml'; path.write_text(yaml.safe_dump(BASE))
    config=load_config(path)
    assert config['run']['sample_id']=='S1'; assert Path(config['inputs']['assembly_fasta']).is_absolute()


def test_paths_resolve_from_configuration_directory(tmp_path, monkeypatch):
    config_dir = tmp_path / "configuration"
    config_dir.mkdir()
    relative_paths = {
        "assembly_fasta": "inputs/assembly.fa",
        "mini_bam": "inputs/sample.bam",
        "mini_bam_index": "inputs/sample.bam.bai",
        "reference_fasta": "reference/reference.fa",
    }
    for relative in relative_paths.values():
        destination = config_dir / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text("fixture\n", encoding="utf-8")

    absolute_index = tmp_path / "absolute-reference.fa.fai"
    absolute_index.write_text("fixture\n", encoding="utf-8")
    locus = config_dir / "loci/locus.yaml"
    locus.parent.mkdir()
    locus.write_text("locus: fixture\n", encoding="utf-8")

    document = deepcopy(BASE)
    document["inputs"].update(relative_paths)
    document["inputs"]["reference_fasta_index"] = str(absolute_index)
    document["locus_config"] = "loci/locus.yaml"
    document["run"]["output_root"] = "results"
    config_path = config_dir / "run.yaml"
    config_path.write_text(yaml.safe_dump(document), encoding="utf-8")

    elsewhere = tmp_path / "unrelated-working-directory"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)
    config = load_config(config_path, check_inputs=True)

    for key, relative in relative_paths.items():
        assert config["inputs"][key] == str((config_dir / relative).resolve())
    assert config["inputs"]["reference_fasta_index"] == str(absolute_index)
    assert config["locus_config"] == str(locus.resolve())
    assert config["run"]["output_root"] == str((config_dir / "results").resolve())

    manifest = create_manifest(config, config_path)
    assert manifest["resolved_paths"] == {
        **config["inputs"],
        "locus_config": config["locus_config"],
        "output_root": config["run"]["output_root"],
    }

def test_reject_missing_required_field(tmp_path):
    broken={**BASE,"run":{"sample_id":"S1","output_root":"outputs"}}; path=tmp_path/'bad.yaml'; path.write_text(yaml.safe_dump(broken))
    with pytest.raises(ConfigurationError, match='run.locus_id'): load_config(path)

@pytest.mark.parametrize('value',['sample_1','HTT-1','A.b'])
def test_safe_identifiers(value): assert validate_identifier(value)==value
@pytest.mark.parametrize('value',['../bad','has space',''])
def test_unsafe_identifiers(value):
    with pytest.raises(ConfigurationError): validate_identifier(value)
