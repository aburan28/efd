#!/usr/bin/env python3
"""Reject stale indices, mismatched evidence and malformed arithmetic programs."""
import ast
import json
from pathlib import Path
try:
    from .build_database import ROOT,build
    from .efd_format import normalize,compile_op3
except ImportError:
    from build_database import ROOT,build
    from efd_format import normalize,compile_op3


def validate():
    path=ROOT/'database/index.json';saved=json.loads(path.read_text());current=build()
    if saved!=current:raise ValueError('Index is stale; run tools/build_database.py')
    models={r['id'] for r in saved['models']};coords={r['id'] for r in saved['coordinates']}
    if len(models)!=len(saved['models']) or len(coords)!=len(saved['coordinates']):raise ValueError('Duplicate model or coordinate IDs')
    ids=set()
    for rec in saved['formulas']:
        if rec['id'] in ids:raise ValueError('Duplicate formula ID')
        ids.add(rec['id'])
        if rec['model'] not in models or rec['coordinates'] not in coords:raise ValueError('Broken coordinate reference')
        if rec['field']=='binary' and not rec['model'].startswith('binary/'):raise ValueError('Field namespace mismatch')
        if not rec['sources']:raise ValueError('Missing source '+rec['id'])
        text=normalize((ROOT/rec['path']).read_text())
        compile_op3([l[8:] for l in text.splitlines() if l.startswith('compute ')])
        for line in (ROOT/rec['op3_path']).read_text().splitlines():
            if not line.strip():continue
            tree=ast.parse(line.replace('^','**')).body
            if len(tree)!=1 or not isinstance(tree[0],ast.Assign):raise ValueError('Malformed operation '+rec['id'])
            for node in ast.walk(tree[0].value):
                if not isinstance(node,(ast.Name,ast.Load,ast.Constant,ast.BinOp,ast.UnaryOp,ast.Add,ast.Sub,ast.Mult,ast.Div,ast.Pow,ast.USub,ast.UAdd)):
                    raise ValueError('Non-arithmetic operation '+rec['id'])
    return saved['counts']


if __name__=='__main__':print(json.dumps(validate(),indent=2))
