#!/usr/bin/env python3
"""Check imported metadata and compiled op3 against independent Weierstrass arithmetic."""
import argparse
from collections import defaultdict
from functools import lru_cache
import json
import hashlib
from pathlib import Path
import re
try:
    from .efd_format import evaluate
except ImportError:
    from efd_format import evaluate

ROOT=Path(__file__).resolve().parents[1]
PRIMES=(101,103,109,127,137,157)


def assignments(path,tag):
    return [(l[len(tag)+1:].split('=',1)[0].strip(),l.split('=',1)[1].lstrip('= ').strip())
            for l in path.read_text().splitlines() if l.startswith(tag+' ')]


def names(s):return set(re.findall(r'[A-Za-z_][A-Za-z_0-9]*',s))


def solve(constraints,env,p):
    env=dict(env);todo=list(constraints)
    for _ in range(len(todo)+2):
        remaining=[];progress=False
        for lhs,rhs in todo:
            unknown=(names(lhs)|names(rhs))-env.keys()
            if len(unknown)==1:
                key=next(iter(unknown))
                roots=[]
                for v in range(p):
                    try:
                        if evaluate(lhs,dict(env,**{key:v}),p)==evaluate(rhs,dict(env,**{key:v}),p):roots.append(v)
                    except (ValueError,ZeroDivisionError):pass
                if not roots:return None
                env[key]=next((v for v in roots if v),roots[0]);progress=True
            elif not unknown:
                try:
                    if evaluate(lhs,env,p)!=evaluate(rhs,env,p):return None
                except (ValueError,ZeroDivisionError):return None
            else:remaining.append((lhs,rhs))
        todo=remaining
        if not todo:return env
        if not progress:return None
    return None


def add(P,Q,a,p):
    if P is None:return Q
    if Q is None:return P
    x,y=P;u,v=Q
    if x==u and (y+v)%p==0:return None
    slope=((3*x*x+a)*pow(2*y,-1,p) if P==Q else (v-y)*pow(u-x,-1,p))%p
    X=(slope*slope-x-u)%p
    return X,(slope*(x-X)-y)%p


@lru_cache(None)
def model_data(model):
    path=ROOT/model/'coordinates'
    params=[l.split()[1] for l in path.read_text().splitlines() if l.startswith('parameter ')]
    coords=[l.split()[1] for l in path.read_text().splitlines() if l.startswith('coordinate ')]
    coeff={}
    for line in path.read_text().splitlines():
        if re.match('a[013246] =',line):
            k,v=line.split('=',1);coeff[k.strip()]=v.strip()
    return params,coords,assignments(path,'satisfying'),assignments(path,'toweierstrass'),coeff


def weier_coeff(model,env,p):
    coeff={k:evaluate(v,env,p) for k,v in model_data(model)[4].items()}
    a0,a1,a2,a3,a4,a6=(coeff[k] for k in ('a0','a1','a2','a3','a4','a6'))
    h2=(a0*a2+a1*a1*pow(4,-1,p))%p
    h4=(a0*a4+a1*a3*pow(2,-1,p))%p
    h6=(a0*a6+a3*a3*pow(4,-1,p))%p
    aa=(a0*h4-h2*h2*pow(3,-1,p))%p
    bb=(2*h2**3*pow(27,-1,p)-h2*a0*h4*pow(3,-1,p)+a0*a0*h6)%p
    if (4*aa**3+27*bb*bb)%p==0:raise ValueError('singular')
    return coeff,h2,aa,bb


def to_weier(model,affine,params,p):
    neutral={'jquartic':{'x':0,'y':1},'jintersect':{'s':0,'c':1,'d':1},
             'twistedhessian':{'x':0,'y':p-1}}
    if model in neutral and affine==neutral[model]:return None
    coeff,h2,aa,bb=weier_coeff(model,params,p)
    env=dict(params,**affine);transformed={}
    for lhs,rhs in model_data(model)[3]:transformed[lhs]=evaluate(rhs,env,p)
    vals=list(transformed.values());u,v=vals
    a0,a1,a3=(coeff[k] for k in ('a0','a1','a3'))
    X=(a0*u+h2*pow(3,-1,p))%p
    Y=a0*(a0*v+(a1*u+a3)*pow(2,-1,p))%p
    if (Y*Y-X**3-aa*X-bb)%p:raise AssertionError(('model map not on curve',model,affine,params,p))
    return X,Y


@lru_cache(None)
def point_pool(model,parameters,p):
    env=dict(parameters);coords=model_data(model)[1];relations=model_data(model)[2]
    points=[]
    if len(coords)==3:
        square=defaultdict(list)
        for t in range(p):square[t*t%p].append(t)
        for s in range(p):
            for c in square[(1-s*s)%p]:
                for d in square[(1-env['a']*s*s)%p]:points.append(dict(s=s,c=c,d=d))
    else:
        for x in range(p):
            for y in range(p):
                e=dict(env,x=x,y=y)
                if all(evaluate(l,e,p)==evaluate(r,e,p) for l,r in relations):points.append(dict(x=x,y=y))
    good=[]
    for point in points:
        try:
            W=to_weier(model,point,env,p)
            if W is not None:good.append((point,W))
        except (ValueError,ZeroDivisionError):pass
    return good


def lift(model,coord,affine,index,params,p,unit=None):
    weights=2 if model in ('3dik','shortw') and coord.startswith('jacobian') else 1
    if model=='3dik':weights=2
    scale=1 if unit=='Z' else (3+2*index)
    if unit and unit!='Z':
        af=unit.lower()
        scale=pow(affine[af],-1,p)
    if model=='jintersect': e={v.upper():affine[v]*scale%p for v in ('s','c','d')}
    else:
        yw=3 if weights==2 else (2 if model in ('2dik','jquartic') else 1)
        e=dict(X=affine['x']*pow(scale,weights,p)%p,Y=affine['y']*pow(scale,yw,p)%p)
    e['Z']=scale
    relations=assignments(ROOT/model/coord/'variables','satisfying')
    todo=[(l,r) for l,r in relations if l not in model_data(model)[1]]
    for _ in range(len(todo)+1):
        nexttodo=[]
        for l,r in todo:
            try:e[l]=evaluate(r,dict(params,**e),p)
            except KeyError:nexttodo.append((l,r))
        todo=nexttodo
    return {k+str(index):v for k,v in e.items()}


def decode(model,coord,result,params,p,index=3):
    variables=[l.split()[1] for l in (ROOT/model/coord/'variables').read_text().splitlines() if l.startswith('variable ')]
    e={k:result[k+str(index)] for k in variables if k+str(index) in result}
    if e.get('Z')==0:raise ZeroDivisionError('output outside affine chart')
    relations=assignments(ROOT/model/coord/'variables','satisfying')
    for _ in range(len(relations)+1):
        for lhs,rhs in relations:
            try:
                v=evaluate(rhs,dict(params,**e),p)
                if lhs in e and lhs not in model_data(model)[1] and e[lhs]!=v:raise AssertionError(('cached-coordinate invariant',lhs,e[lhs],v))
                e[lhs]=v
            except KeyError:pass
    affine={k:e[k] for k in model_data(model)[1]}
    if any(evaluate(l,dict(params,**affine),p)!=evaluate(r,dict(params,**affine),p) for l,r in model_data(model)[2]):raise AssertionError(('output off model',affine))
    return to_weier(model,affine,params,p)


def execute(lines,env,p):
    env=dict(env)
    for lhs,rhs in lines:env[lhs]=evaluate(rhs,env,p)
    return env


def verify_formula(path):
    model,coord,operation,_=path.split('/');full=ROOT/path
    formula_assume=assignments(full,'assume')
    coord_assume=assignments(ROOT/model/coord/'variables','assume')
    units={int(l[-1]):l[:-1] for l,r in formula_assume if re.fullmatch('[A-Z][12]',l) and r=='1'}
    constraints=coord_assume+[(l,r) for l,r in formula_assume if not re.search(r'\b[A-Z][A-Z0-9]*[12]\b',l+r)]
    params=model_data(model)[0];checks=0;fields=set();configs=0;exceptional=0
    seen_configs=set();seen_cases=set()
    programs=[assignments(full,'compute')]
    programs.append([(l.split('=',1)[0].strip(),l.split('=',1)[1].strip()) for l in full.with_suffix('.op3').read_text().splitlines()])
    for p in PRIMES:
        for av in (2,3,5,7,11):
            env={k:(av if k in ('a','d') else 2) for k in params}
            if model=='shortw':env['b']=5
            for lhs,rhs in coord_assume:
                if lhs in params and names(rhs)<=env.keys():env[lhs]=evaluate(rhs,env,p)
            env=solve(constraints,env,p)
            if env is None:continue
            config_key=(p,tuple(sorted(env.items())))
            if config_key in seen_configs:continue
            seen_configs.add(config_key)
            base={k:env[k] for k in params}
            try:
                _,_,aa,_=weier_coeff(model,base,p)
                pool=point_pool(model,tuple(sorted(base.items())),p)
            except (ValueError,ZeroDivisionError):continue
            if len(pool)<8:continue
            perconfig=0
            for i in range(12):
                P,WP=pool[(7*i+3)%len(pool)];Q,WQ=pool[(13*i+7)%len(pool)]
                try:
                    inp=dict(env,**lift(model,coord,P,1,env,p,units.get(1)))
                    if operation=='addition':inp.update(lift(model,coord,Q,2,env,p,units.get(2)))
                    if solve(formula_assume,inp,p) is None:continue
                    case_key=(p,tuple(sorted(inp.items())))
                    if case_key in seen_cases:continue
                    seen_cases.add(case_key)
                    expected=(add(WP,WQ,aa,p) if operation=='addition' else
                              add(add(WP,WP,aa,p),WP,aa,p) if operation=='tripling' else
                              add(WP,WP,aa,p) if operation=='doubling' else WP)
                    results=[execute(pr,inp,p) for pr in programs]
                    # Independent execution of metadata and lowered operation code.
                    for result in results:
                        try:actual=decode(model,coord,result,env,p)
                        except ZeroDivisionError:raise ZeroDivisionError('output chart excluded')
                        if actual!=expected:raise AssertionError(('group mismatch',p,env,P,Q,actual,expected))
                    outputs=set(k for k in results[0] if re.search('3$',k))
                    if any(results[0][k]!=results[1].get(k) for k in outputs):raise AssertionError('op3 mismatch')
                except (ZeroDivisionError,ValueError):exceptional+=1;continue
                checks+=1;perconfig+=1;fields.add(p)
            if perconfig:configs+=1
            if checks>=24 and len(fields)>=2:return dict(status='finite_field_checked',cases=checks,fields=sorted(fields),parameter_sets=configs,chart_or_denominator_exclusions=exceptional)
    return dict(status='insufficient_coverage',cases=checks,fields=sorted(fields),parameter_sets=configs,chart_or_denominator_exclusions=exceptional)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'research/import-verification.json')
    args=parser.parse_args()
    manifest=json.loads((ROOT/'research/import-manifest.json').read_text())
    records=[]
    for rec in manifest['files']:
        if not rec['operation_code']:continue
        try:result=verify_formula(rec['path'])
        except Exception as e:result=dict(status='failed',error=repr(e))
        records.append(dict(path=rec['path'],metadata_sha256=hashlib.sha256((ROOT/rec['path']).read_bytes()).hexdigest(),op3_sha256=hashlib.sha256((ROOT/(rec['path']+'.op3')).read_bytes()).hexdigest(),**result))
        if result['status']!='finite_field_checked':print(rec['path'],result,flush=True)
    report=dict(method='Exact arithmetic over multiple small prime fields; metadata and compiled op3 independently executed and compared with a generalized-Weierstrass-to-short-Weierstrass group law. On affine charts only; no completeness or symbolic proof claimed.',formulas=records)
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    counts=defaultdict(int)
    for rec in records:counts[rec['status']]+=1
    print(dict(counts),flush=True)
    if counts.get('failed') or counts.get('insufficient_coverage'):raise SystemExit(1)


if __name__=='__main__':main()
