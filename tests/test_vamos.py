from pathlib import Path
import pytest
from tr_calling_pipeline.callers.vamos import UnsupportedVamosFormat, UnsupportedVamosVersion, parse_native_jsonl, select_adapter

FIX=Path(__file__).parent/"fixtures"/"vamos"
def test_supported_adapter_and_unknown_rejected():
    assert select_adapter("2.1.0").capabilities.supports_bam_input
    with pytest.raises(UnsupportedVamosVersion): select_adapter("3.0.0")
    with pytest.raises(UnsupportedVamosVersion): select_adapter(None)
def test_lossless_known_parser():
    records=parse_native_jsonl(FIX/"supported-version"/"read.native.jsonl")
    assert [r["allele_id"] for r in records]==["allele-1","allele-2"]
    assert records[0]["repeat_count"]=="19"
    assert records[0]["motif_chain"][1]["motif"]=="CAA"
def test_malformed_rejected():
    with pytest.raises(UnsupportedVamosFormat): parse_native_jsonl(FIX/"malformed"/"malformed.jsonl")
