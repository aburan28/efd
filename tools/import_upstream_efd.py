#!/usr/bin/env python3
"""Stage missing prime-field EFD data, normalize products, and compile clean op3 code."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import date
import hashlib
import json
from pathlib import Path
import re
from urllib.request import urlopen
from efd_format import normalize, compile_op3

ROOT=Path(__file__).resolve().parents[1]
BASE='https://www.hyperelliptic.org/EFD/g1p/'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache',type=Path,required=True)
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    snapshot=json.loads((ROOT/'research/upstream-efd-inventory.json').read_text())
    paths=snapshot['missing_models']+snapshot['missing_coordinates']+snapshot['missing_formulas']
    def fetch(path):
        cached=args.cache/path
        if cached.exists():raw=cached.read_bytes()
        else:
            with urlopen(BASE+'data/'+path,timeout=30) as r:raw=r.read()
            cached.parent.mkdir(parents=True,exist_ok=True);cached.write_bytes(raw)
        normalized=normalize(raw.decode())
        if path.endswith('/coordinates'):
            normalized='\n'.join(l.replace(' = ',' == ') if l.startswith('satisfying ') and ' == ' not in l else l for l in normalized.splitlines())+'\n'
        if '/diffadd/' in path:
            normalized='\n'.join(re.sub(r'\b([A-Z][A-Z]*)5\b',r'\g<1>4',l) if l.startswith('compute ') else l for l in normalized.splitlines())+'\n'
        computations=[l[8:] for l in normalized.splitlines() if l.startswith('compute ')]
        op3=compile_op3(computations) if computations else None
        return path,raw,normalized,op3
    manifest=[]
    with ThreadPoolExecutor(max_workers=4) as pool:
        for path,raw,normalized,op3 in pool.map(fetch,paths):
            record=dict(path=path,url=BASE+'data/'+path,sha256=hashlib.sha256(raw).hexdigest(),
                        normalized_sha256=hashlib.sha256(normalized.encode()).hexdigest(),
                        verification='pending',operation_code='compiled_from_metadata' if op3 else None)
            if args.write:
                dest=ROOT/path
                if dest.exists():
                    existing=dest.read_text()
                    enriched=path.endswith(('/coordinates','/variables')) and existing.startswith(normalized)
                    if existing!=normalized and not enriched:raise ValueError('Refusing to overwrite changed file '+path)
                else:
                    dest.parent.mkdir(parents=True,exist_ok=True);dest.write_text(normalized)
                if op3:
                    code=dest.with_suffix('.op3')
                    if code.exists() and code.read_text()!=op3:raise ValueError('Refusing to overwrite changed operation code '+path)
                    code.write_text(op3)
            manifest.append(record)
    if args.write:
        (ROOT/'research/import-manifest.json').write_text(json.dumps(dict(imported_on=date.today().isoformat(),source=BASE,
            policy='Implicit products normalized; operation code regenerated from metadata, preserving source expressions. Verification recorded by a separate verifier.',files=manifest),indent=2)+'\n')
    print(json.dumps(dict(data_files=len(manifest),formula_files=sum(r['operation_code'] is not None for r in manifest),written=args.write)))


if __name__=='__main__':main()
