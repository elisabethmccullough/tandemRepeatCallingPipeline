#!/usr/bin/env python3
"""Cross-platform synthetic executable for the provisional adapter."""
import os, pathlib, sys
version=os.environ.get("FAKE_STRAGLR_VERSION","STRaglr 1.5.0")
if "--version" in sys.argv:
    target=sys.stderr if os.environ.get("FAKE_STRAGLR_VERSION_STDERR") else sys.stdout
    print(version,file=target); raise SystemExit(0)
if os.environ.get("FAKE_STRAGLR_FAIL"):
    print("synthetic failure",file=sys.stderr); raise SystemExit(7)
if len(sys.argv)<5: raise SystemExit(2)
bam,reference,catalog,output=sys.argv[1:5]
if os.environ.get("FAKE_STRAGLR_ECHO"): print("\n".join(sys.argv[1:]))
if os.environ.get("FAKE_STRAGLR_OMIT"): raise SystemExit(0)
p=pathlib.Path(output); p.parent.mkdir(parents=True,exist_ok=True)
if os.environ.get("FAKE_STRAGLR_EMPTY"): p.write_text(""); raise SystemExit(0)
if os.environ.get("FAKE_STRAGLR_UNSUPPORTED"): p.write_text("unknown\tcolumns\nx\ty\n"); raise SystemExit(0)
if os.environ.get("FAKE_STRAGLR_MALFORMED"): p.write_text("locus_id\tchromosome\tstart\tend\tallele_id\nL\tchr1\t1\t2\n"); raise SystemExit(0)
p.write_text("record_id\tlocus_id\tchromosome\tstart\tend\tallele_id\trepeat_unit\trepeat_count\trepeat_size\tsupporting_reads\ttotal_spanning_reads\tcaller_specific\nrecord-1\tSYNTH_LOCUS\tchrSynthetic\t100\t112\tallele-1\tCAG\t004\t12\t3\t8\t0007\nrecord-2\tSYNTH_LOCUS\tchrSynthetic\t100\t118\tallele-2\tCAG\t6\t18\t5\t8\talpha\n")
if os.environ.get("FAKE_STRAGLR_MULTIPLE"):
    p.with_suffix(".details.txt").write_text("synthetic auxiliary native output\n")
