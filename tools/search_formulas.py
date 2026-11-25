#!/usr/bin/env python3
"""Search native formulas without executing their expressions."""
import argparse
from collections import Counter
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--field',choices=['prime','binary'])
    p.add_argument('--model',help='Model ID or unnamespaced model name')
    p.add_argument('--operation')
    p.add_argument('--source',help='Case-insensitive source substring')
    p.add_argument('--verified',action='store_true',help='Only independently finite-field-checked entries')
    p.add_argument('--json',action='store_true')
    p.add_argument('--summary',action='store_true')
    a=p.parse_args();db=json.loads((ROOT/'database/index.json').read_text())
    rows=[f for f in db['formulas'] if (not a.field or f['field']==a.field)
          and (not a.model or f['model']==a.model or f['model'].split('/')[-1]==a.model)
          and (not a.operation or f['operation']==a.operation)
          and (not a.source or a.source.lower() in ' '.join(f['sources']).lower())
          and (not a.verified or f['verification']['status']=='finite_field_checked')]
    if a.json:print(json.dumps(rows,indent=2))
    elif a.summary:print(json.dumps(dict(entries=len(rows),by_field=dict(Counter(f['field'] for f in rows)),by_model=dict(Counter(f['model'] for f in rows)),by_operation=dict(Counter(f['operation'] for f in rows)),verification=dict(Counter(f['verification']['status'] for f in rows))),indent=2))
    else:
        for f in rows:print(f['id']+' ['+f['verification']['status']+']')
        print(str(len(rows))+' entries')


if __name__=='__main__':main()
