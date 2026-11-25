#!/usr/bin/env python3
"""Build the unified, field-aware index from native formula files and verification reports."""
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def tags(path):
    data={}
    for line in path.read_text().splitlines():
        key,sep,value=line.partition(' ')
        if sep:data.setdefault(key,[]).append(value)
        elif line.strip():data.setdefault(key,[])
    return data


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def build():
    atlas=json.loads((ROOT/'transformations/catalog.json').read_text())
    verification={}
    for name in ('import-verification.json','binary-verification.json','additional-formulas-verification.json'):
        p=ROOT/'research'/name
        if p.exists():
            for rec in json.loads(p.read_text())['formulas']:
                verification[rec['path']]=dict(rec,report='research/'+name)
    paper={rec['path']:rec for rec in json.loads((ROOT/'research/paper-imports.json').read_text())['formulas']}
    if any(r['verification']=='pending' for r in paper.values()):raise ValueError('Paper verification status must be recorded before indexing')
    manifest={}
    for name in ('import-manifest.json','binary-import-manifest.json'):
        for rec in json.loads((ROOT/'research'/name).read_text())['files']:manifest[rec['path']]=rec
    formulas=[];models=[];coords=[]
    model_paths=sorted(ROOT.glob('*/coordinates'))+sorted((ROOT/'binary').glob('*/coordinates'))
    for path in model_paths:
        model=path.parent;model_id=str(model.relative_to(ROOT));field='binary' if model_id.startswith('binary/') else 'prime'
        m=tags(path)
        # Native a-invariants have lines like a0 = 1; retain their exact expressions.
        invariants={k:v[0].lstrip('= ') for k,v in m.items() if k in ('a0','a1','a2','a3','a4','a6')}
        nonsingular=atlas['models'].get(model_id,{}).get('nonsingular')
        if field=='binary':nonsingular={'shortw':'a6 != 0','edwards':'d1 != 0 and d2 != d1^2+d1','hessian':'d != 0 and d^3 != 1 (ordinary binary family)'}[model.name]
        models.append(dict(id=model_id,field=field,characteristic=2 if field=='binary' else 'greater_than_3',path=str(path.relative_to(ROOT)),name=m['name'][0],parameters=m.get('parameter',[]),affine_coordinates=m.get('coordinate',[]),equations=m.get('satisfying',[]),nonsingular=nonsingular,weierstrass_map=dict(forward=m.get('toweierstrass',[]),inverse=m.get('fromweierstrass',[]),invariants=invariants),source=manifest.get(str(path.relative_to(ROOT)),{}).get('url')))
        for var in sorted(model.glob('*/variables')):
            coord=var.parent;c=tags(var);coord_id=str(coord.relative_to(ROOT))
            coords.append(dict(id=coord_id,model=model_id,field=field,path=str(var.relative_to(ROOT)),name=c['name'][0],variables=c.get('variable',[]),parameters=c.get('parameter',[]),assumptions=c.get('assume',[]),relations=c.get('satisfying',[]),to_affine=c.get('toaffine',[]),from_affine=c.get('tosystem',[])))
            for op3 in sorted(coord.glob('*/*.op3')):
                meta=op3.with_suffix('');rel=str(meta.relative_to(ROOT));f=tags(meta)
                v=verification.get(rel)
                if v:
                    if v.get('metadata_sha256')!=digest(meta) or v.get('op3_sha256')!=digest(op3):
                        raise ValueError('Stale verification evidence for '+rel)
                    status=v['status'];evidence={k:v[k] for k in ('report','cases','fields') if k in v}
                elif rel in paper:
                    if paper[rel].get('metadata_sha256')!=digest(meta) or paper[rel].get('op3_sha256')!=digest(op3):
                        raise ValueError('Stale paper verification evidence for '+rel)
                    status=paper[rel]['verification'];evidence={'test':'tests/test_paper_imports.py'}
                else:status='inherited_not_rechecked';evidence={}
                sources=f.get('source',[]) or [manifest.get(rel,{}).get('url') or 'https://github.com/aburan28/efd/blob/e82e4dfddbcf866f18992beeed622e5c0083ca63/'+rel]
                formulas.append(dict(id=rel,path=rel,op3_path=str(op3.relative_to(ROOT)),model=model_id,coordinates=coord_id,field=field,operation=meta.parent.name,sources=sources,attribution='paper_source' if f.get('source') else 'metadata_or_export_source_only',parameters=f.get('parameter',[]),assumptions=c.get('assume',[])+f.get('assume',[]),unified='unified' in f,verification=dict(status=status,**evidence),sha256=digest(meta),op3_sha256=digest(op3),contract=paper.get(rel,{}).get('contract')))
    counts=dict(model_families=len(models),coordinate_representations=len(coords),formula_entries=len(formulas),prime_model_families=sum(m['field']=='prime' for m in models),binary_model_families=sum(m['field']=='binary' for m in models),prime_formulas=sum(f['field']=='prime' for f in formulas),binary_formulas=sum(f['field']=='binary' for f in formulas))
    return dict(schema_version=1,scope='Genus-one formulas with explicit prime/binary namespaces, model equations and native Weierstrass maps. Counts are entries, not algebraically independent identities.',domain_policy='Equality assumptions are retained; nonzero-denominator and completeness conditions are not universally encoded. Consult the source and explicit operation contracts. Unified never implies complete.',verification_semantics={'finite_field_checked':'Checked on stated finite fields and affine domains, not a proof of completeness or symbolic correctness for every field.','inherited_not_rechecked':'Existing export entry retained; not covered by the new independent arithmetic checks.','syntax_checked_only':'Parsed and compiled, not yet checked against group arithmetic.'},counts=counts,models=models,coordinates=coords,formulas=formulas)


if __name__=='__main__':
    result=build();dest=ROOT/'database/index.json';dest.parent.mkdir(exist_ok=True);dest.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result['counts'],indent=2))
