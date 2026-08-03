#!/usr/bin/env python3
"""Portable VAMOS protocol double for adapter tests; not a biological caller."""
import json, os, pathlib, sys
args=sys.argv[1:]
if "--version" in args:
    stream=sys.stderr if os.environ.get("FAKE_VAMOS_VERSION_STDERR") else sys.stdout
    print(os.environ.get("FAKE_VAMOS_VERSION","VAMOS 2.1.0"),file=stream); raise SystemExit(0)
if os.environ.get("FAKE_VAMOS_FAIL"): raise SystemExit(17)
prefix=pathlib.Path(args[args.index("--output-prefix")+1]); prefix.parent.mkdir(parents=True,exist_ok=True)
mode="read" if "--bam" in args else "contig"
input_path=pathlib.Path(args[args.index("--bam" if mode=="read" else "--fasta")+1])
if os.environ.get("FAKE_VAMOS_UNSUPPORTED"):
    prefix.with_suffix(".jsonl").write_text("unsupported\n")
else:
    identifiers=["allele-1","allele-2"] if mode=="read" else ["contig-consensus"]
    header=input_path.read_text().splitlines()[0][1:] if mode=="contig" else None
    with prefix.with_suffix(".jsonl").open("w") as out:
        for i, allele in enumerate(identifiers,1):
            json.dump({"record_id":"synthetic-locus","allele_id":allele,"reference_build":"GRCh38","chromosome":"chr4","start":10,"end":20,"coordinate_convention":"zero_based_half_open","motif":"CAG","motif_chain":[{"motif":"CAG","count":str(10*i)}],"repeat_count":str(10*i),"repeat_length_bp":str(30*i),"supporting_reads":i if mode=="read" else None,"total_spanning_reads":3 if mode=="read" else None,"quality_state":"AVAILABLE" if mode=="read" else "NOT_APPLICABLE","source_header":header},out); out.write("\n")
if not os.environ.get("FAKE_VAMOS_OMIT_SUMMARY"): prefix.with_suffix(".summary.txt").write_text(f"mode={mode}\ninput={input_path}\n")
