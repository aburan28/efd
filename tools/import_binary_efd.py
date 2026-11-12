#!/usr/bin/env python3
"""Import the complete ordinary-binary EFD into an explicit characteristic-2 namespace."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import date
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import re
from urllib.request import urlopen
from efd_format import normalize,compile_op3

ROOT=Path(__file__).resolve().parents[1]
BASE='https://www.hyperelliptic.org/EFD/g12o/'


class Links(HTMLParser):
    def __init__(self,s):super().__init__();self.links=set();self.feed(s)
    def handle_starttag(self,t,a):
        if t=='a':self.links.add(dict(a).get('href',''))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--cache',type=Path,required=True);p.add_argument('--write',action='store_true');args=p.parse_args()
    def get(name):
        cached=args.cache/name
        if cached.exists():return cached.read_bytes()
        with urlopen(BASE+name,timeout=30) as resp:data=resp.read()
        cached.parent.mkdir(parents=True,exist_ok=True);cached.write_bytes(data);return data
    pages=sorted(h for h in Links(get('index.html').decode()).links if h.startswith('auto-') and h.endswith('.html'))
    refs=set()
    with ThreadPoolExecutor(max_workers=4) as pool:
        for raw in pool.map(get,pages):refs.update(h[5:] for h in Links(raw.decode()).links if h.startswith('data/'))
    records=[]
    def stage(path):
        raw=get('data/'+path);text=normalize(raw.decode())
        if path.endswith('/coordinates'):text='\n'.join(l.replace(' = ',' == ') if l.startswith('satisfying ') and ' == ' not in l else l for l in text.splitlines())+'\n'
        if '/diffadd/' in path:text='\n'.join(re.sub(r'\b([A-Za-z][A-Za-z]*)5\b',r'\g<1>4',l) if l.startswith('compute ') else l for l in text.splitlines())+'\n'
        computations=[l[8:] for l in text.splitlines() if l.startswith('compute ')]
        op3=compile_op3(computations) if computations else None
        dest=ROOT/'binary'/path
        if args.write:
            if dest.exists() and dest.read_text()!=text:raise ValueError('Refusing to overwrite changed binary file '+path)
            dest.parent.mkdir(parents=True,exist_ok=True);dest.write_text(text)
            if op3:dest.with_suffix('.op3').write_text(op3)
        return dict(path='binary/'+path,url=BASE+'data/'+path,sha256=hashlib.sha256(raw).hexdigest(),normalized_sha256=hashlib.sha256(text.encode()).hexdigest(),verification='syntax_checked_only',operation_code='compiled_from_metadata' if op3 else None)
    with ThreadPoolExecutor(max_workers=4) as pool:
        for record in pool.map(stage,sorted(refs)):records.append(record)
    result=dict(imported_on=date.today().isoformat(),source=BASE,characteristic=2,
                conditions='Binary field formulas only. Do not evaluate using arithmetic modulo an odd prime. Trace/root assumptions remain necessary.',files=records)
    if args.write:(ROOT/'research/binary-import-manifest.json').write_text(json.dumps(result,indent=2)+'\n')
    print(dict(models=sum(r['path'].endswith('/coordinates') for r in records),representations=sum(r['path'].endswith('/variables') for r in records),formulas=sum(bool(r['operation_code']) for r in records),written=args.write))


if __name__=='__main__':main()
