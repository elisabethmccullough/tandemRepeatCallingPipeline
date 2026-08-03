#!/usr/bin/env python3
import os,sys
from pathlib import Path
name=Path(sys.argv[0]).name
if '--version' in sys.argv:
 print(os.getenv('FAKE_VERSION') or ('LAST 1450' if name in ('lastdb','lastal') else 'tandem-genotypes 0.1.0'),file=sys.stderr if os.getenv('VERSION_STDERR') else sys.stdout); raise SystemExit(0)
if os.getenv('FAKE_EXIT'): raise SystemExit(int(os.getenv('FAKE_EXIT')))
if name=='lastdb':
 prefix=Path(sys.argv[-2]); fasta=Path(sys.argv[-1]); prefix.with_suffix('.prj').write_text('fake database\n'); prefix.with_suffix('.source').write_text(fasta.read_text().splitlines()[0][1:]+'\n')
elif name=='lastal':
 prefix=Path(sys.argv[-2]); seq=prefix.with_suffix('.source').read_text().strip()
 if os.getenv('FAKE_UNSUPPORTED'): print('unsupported')
 else: print(f'##maf version=1\na score=10\ns chr4 0 3 + 10 CAG\ns {seq} 0 3 + 3 CAG')
else:
 out=Path(sys.argv[3]); out.parent.mkdir(parents=True,exist_ok=True)
 if not os.getenv('FAKE_MISSING'):
  out.write_text('' if os.getenv('FAKE_EMPTY') else f'record_id\tallele_id\tsequence_id\tmotif\trepeat_count\nrec-{out.stem}\tallele-1\t{out.stem}\tCAG\t020\n')
