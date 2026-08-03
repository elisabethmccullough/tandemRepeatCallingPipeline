"""Cross-platform protocol doubles used only by the bundled synthetic demo."""
from __future__ import annotations

import json
from pathlib import Path
import shutil
import sys


def _name() -> str:
    return Path(sys.argv[0]).name.lower().replace("_", "-")


def samtools() -> int:
    args = sys.argv[1:]
    if not args or "--version" in args:
        print("samtools 1.20 synthetic fake tool"); return 0
    command = args[0]
    if command == "quickcheck": return 0
    if command == "view" and "-H" in args:
        print("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chrSynthetic\tLN:96"); return 0
    if command == "view":
        output = Path(args[args.index("-o") + 1]); shutil.copyfile(Path(args[-1]), output); return 0
    if command == "sort":
        output = Path(args[args.index("-o") + 1]); shutil.copyfile(Path(args[-1]), output); return 0
    if command == "index":
        output = Path(args[args.index("-o") + 1]); output.write_text("SYNTHETIC_BAI\n"); return 0
    if command == "flagstat":
        print("2 + 0 in total\n0 + 0 secondary\n0 + 0 supplementary\n2 + 0 mapped\n0 + 0 duplicates"); return 0
    if command == "idxstats":
        print("chrSynthetic\t96\t2\t0"); return 0
    return 2


def minimap2() -> int:
    if "--version" in sys.argv[1:]: print("2.28"); return 0
    fasta = Path(sys.argv[-1]); records=[]; name=None; sequence=[]
    for line in fasta.read_text().splitlines():
        if line.startswith(">"):
            if name: records.append((name, "".join(sequence)))
            name=line[1:].split()[0]; sequence=[]
        else: sequence.append(line)
    if name: records.append((name, "".join(sequence)))
    print("@HD\tVN:1.6\tSO:unknown\n@SQ\tSN:chrSynthetic\tLN:96")
    for index, (identifier, sequence) in enumerate(records):
        print(f"{identifier}\t0\tchrSynthetic\t{10 + index * 20}\t60\t{len(sequence)}M\t*\t0\t0\t{sequence}\t*")
    return 0


def vamos() -> int:
    args=sys.argv[1:]
    if "--version" in args: print("VAMOS 2.1.0 synthetic fake tool"); return 0
    prefix=Path(args[args.index("--output-prefix")+1]); prefix.parent.mkdir(parents=True,exist_ok=True)
    mode="read" if "--bam" in args else "contig"
    source_header=None
    if mode == "contig": source_header=Path(args[args.index("--fasta")+1]).read_text().splitlines()[0][1:]
    identifiers=("allele-1","allele-2") if mode == "read" else ("contig-call",)
    with prefix.with_suffix(".jsonl").open("w", encoding="utf-8") as stream:
        for index, allele in enumerate(identifiers, 1):
            json.dump({"record_id":source_header or "synthetic-locus","allele_id":allele,"reference_build":"SYNTHETIC",
                "chromosome":"chrSynthetic","start":9,"end":27,"coordinate_convention":"zero_based_half_open",
                "motif":"CAG","motif_chain":[{"motif":"CAG","count":str(index+3)}],
                "repeat_count":str(index+3),"repeat_length_bp":str((index+3)*3),
                "supporting_reads":index if mode=="read" else None,"total_spanning_reads":2 if mode=="read" else None,
                "quality_state":"AVAILABLE" if mode=="read" else "NOT_APPLICABLE","source_header":source_header}, stream)
            stream.write("\n")
    prefix.with_suffix(".summary.txt").write_text(f"synthetic_mode={mode}\n")
    return 0


def straglr() -> int:
    if "--version" in sys.argv[1:]: print("STRaglr 1.5.0 synthetic fake tool"); return 0
    output=Path(sys.argv[4]); output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text("record_id\tlocus_id\tchromosome\tstart\tend\tallele_id\trepeat_unit\trepeat_count\trepeat_size\tsupporting_reads\ttotal_spanning_reads\tcaller_specific\nrecord-1\tSYNTHETIC-LOCUS\tchrSynthetic\t9\t27\tallele-1\tCAG\t4\t12\t2\t2\tsynthetic\n")
    return 0


def lastdb() -> int:
    if "--version" in sys.argv[1:]: print("LAST 1450 synthetic fake tool"); return 0
    prefix=Path(sys.argv[-2]); fasta=Path(sys.argv[-1]); prefix.with_suffix(".prj").write_text("synthetic database\n")
    prefix.with_suffix(".source").write_text(fasta.read_text().splitlines()[0][1:] + "\n"); return 0


def lastal() -> int:
    if "--version" in sys.argv[1:]: print("LAST 1450 synthetic fake tool"); return 0
    prefix=Path(sys.argv[-2]); sequence_id=prefix.with_suffix(".source").read_text().strip()
    print(f"##maf version=1\na score=10\ns chrSynthetic 9 6 + 96 CAGCAG\ns {sequence_id} 0 6 + 6 CAGCAG"); return 0


def tandem_genotypes() -> int:
    if "--version" in sys.argv[1:]: print("tandem-genotypes 0.1.0 synthetic fake tool"); return 0
    alignment, _, output=sys.argv[1:4]; sequence_id=Path(alignment).stem
    Path(output).write_text("record_id\tallele_id\tsequence_id\treference_build\tchromosome\tmotif\trepeat_count\n"
        f"record-{sequence_id}\tallele-1\t{sequence_id}\tSYNTHETIC\tchrSynthetic\tCAG\t4\n"); return 0


def main() -> int:
    name=_name()
    for key, function in (("samtools",samtools),("minimap2",minimap2),("vamos",vamos),("straglr",straglr),
                          ("lastdb",lastdb),("lastal",lastal),("tandem-genotypes",tandem_genotypes)):
        if key in name: return function()
    raise SystemExit(f"unknown synthetic fake tool entry point: {name}")
