import json
from decimal import Decimal
from pathlib import Path
import pytest
import yaml
from tr_calling_pipeline.case_package import validate_case_package
from tr_calling_pipeline.config import ConfigurationError, load_config, load_locus_config
from tr_calling_pipeline.models import AnalysisSource, AssignmentState, CallerEvidence, EvidenceState, SequenceRole

ROOT = Path(__file__).parents[1]

def test_example_run_and_locus_configs_validate():
    run = load_config(ROOT / "config/example.yaml")
    assert Path(run["inputs"]["assembly_fasta"]).is_absolute()
    assert run["inputs"]["assembly_records"][1]["sequence_role"] == "PATIENT_HAPLOTYPE"
    locus = load_locus_config(ROOT / "config/loci/htt_hg38.yaml")
    assert locus["locus"]["target_region"]["start"] is None

def test_invalid_run_identifier(tmp_path):
    data = yaml.safe_load((ROOT / "config/example.yaml").read_text())
    data["run"]["case_id"] = "../escape"
    path = tmp_path / "bad.yaml"; path.write_text(yaml.safe_dump(data))
    with pytest.raises(ConfigurationError): load_config(path)

def test_consensus_role_is_preserved():
    assert SequenceRole.PATIENT_CONSENSUS.value != SequenceRole.PATIENT_HAPLOTYPE.value

def test_duplicate_blocks_and_motifs(tmp_path):
    data = yaml.safe_load((ROOT / "config/loci/htt_hg38.yaml").read_text())
    duplicate = dict(data["repeat_blocks"][0]); data["repeat_blocks"].append(duplicate)
    path = tmp_path / "locus.yaml"; path.write_text(yaml.safe_dump(data))
    with pytest.raises(ConfigurationError, match="block IDs"): load_locus_config(path)
    data["repeat_blocks"] = [data["repeat_blocks"][0]]
    data["repeat_blocks"][0]["motifs"].append(dict(data["repeat_blocks"][0]["motifs"][0]))
    path.write_text(yaml.safe_dump(data))
    with pytest.raises(ConfigurationError, match="motif IDs"): load_locus_config(path)

def test_invalid_coordinate_convention(tmp_path):
    data = yaml.safe_load((ROOT / "config/loci/htt_hg38.yaml").read_text())
    data["measurement"]["coordinate_convention"] = "guess"
    path = tmp_path / "locus.yaml"; path.write_text(yaml.safe_dump(data))
    with pytest.raises(ConfigurationError): load_locus_config(path)

def test_minimal_and_multi_locus_fixtures():
    minimal = validate_case_package(ROOT / "tests/fixtures/case-packages/minimal")
    assert minimal["loci"][0]["patient_sequences"][2]["sequence_role"] == "PATIENT_CONSENSUS"
    assert len(validate_case_package(ROOT / "tests/fixtures/case-packages/multi-locus")["loci"]) == 2

@pytest.mark.parametrize("name", ["invalid-traversal", "invalid-absolute"])
def test_unsafe_path_fixtures_are_invalid(name):
    with pytest.raises(ConfigurationError): validate_case_package(ROOT / "tests/fixtures/case-packages" / name)

def test_duplicate_file_and_sequence_ids(tmp_path):
    source = json.loads((ROOT / "tests/fixtures/case-packages/minimal/case-manifest.json").read_text())
    source["files"].append(dict(source["files"][0]))
    (tmp_path / "case-manifest.json").write_text(json.dumps(source))
    with pytest.raises(ConfigurationError, match="duplicate file_id"): validate_case_package(tmp_path)
    source["files"].pop(); source["loci"][0]["patient_sequences"][1]["sequence_id"] = "patient-hap1"
    (tmp_path / "case-manifest.json").write_text(json.dumps(source))
    with pytest.raises(ConfigurationError, match="duplicate sequence_id"): validate_case_package(tmp_path)

def test_symlink_escape(tmp_path):
    fixture = json.loads((ROOT / "tests/fixtures/case-packages/minimal/case-manifest.json").read_text())
    fixture["files"][0].update(path="escape/secret", required=True, status="AVAILABLE")
    (tmp_path / "case-manifest.json").write_text(json.dumps(fixture))
    outside = tmp_path.parent / "outside"; outside.mkdir(exist_ok=True); (outside / "secret").write_text("x")
    (tmp_path / "escape").symlink_to(outside, target_is_directory=True)
    with pytest.raises(ConfigurationError, match="symlink"): validate_case_package(tmp_path)

def test_unknown_file_reference(tmp_path):
    fixture = json.loads((ROOT / "tests/fixtures/case-packages/minimal/case-manifest.json").read_text())
    fixture["loci"][0]["patient_sequences"][0]["source_fasta_file_id"] = "unknown"
    (tmp_path / "case-manifest.json").write_text(json.dumps(fixture))
    with pytest.raises(ConfigurationError, match="unknown file"): validate_case_package(tmp_path)

def test_lossless_evidence_and_source_separation():
    common = dict(record_schema_version="1.0", record_id="r1", case_id="C", subject_id="S", sample_id="X", locus_id="HTT", caller="VAMOS", caller_version="x", source_file_id="native", source_file_sha256="0"*64, native_record_identifier="line1", native_allele_identifier=None, associated_sequence_id=None, assignment_state=AssignmentState.UNASSIGNED, reference_build="GRCh38", chromosome="chr4", start=None, end=None, coordinate_convention="zero_based_half_open", reported_motif="CAG", reported_motif_chain=({"motif":"CAG","count":"40.50"},), reported_repeat_count=Decimal("40.50"), reported_repeat_length_bp=None, supporting_reads=None, total_spanning_reads=None, quality_state=EvidenceState.AVAILABLE, raw_fields={"decimal":"40.50", "array":["CAG","CAA"]}, normalization_warnings=("length not reported",))
    raw = CallerEvidence(analysis_source=AnalysisSource.RAW_READS, **common)
    assembled = CallerEvidence(analysis_source=AnalysisSource.ASSEMBLED_CONTIG, **{**common,"record_id":"r2","associated_sequence_id":"patient-hap1","assignment_state":AssignmentState.DIRECT_SEQUENCE_ASSOCIATION})
    assert str(raw.reported_repeat_count) == "40.50"
    assert raw.raw_fields["decimal"] == "40.50" and raw.reported_motif_chain[0]["count"] == "40.50"
    assert raw.analysis_source != assembled.analysis_source
