#!/usr/bin/env python3
"""Verify ordinary-binary formulas over polynomial-basis extension fields."""
import argparse
import ast
from collections import Counter
from functools import lru_cache
import json
import hashlib
from pathlib import Path
import re
try:
    from .efd_format import expression_tree
except ImportError:
    from efd_format import expression_tree

ROOT=Path(__file__).resolve().parents[1]


class GF2:
    def __init__(self,modulus):self.modulus=modulus;self.q=1<<(modulus.bit_length()-1)
    def add(self,a,b):return a^b
    def mul(self,a,b):
        r=0
        while b:
            if b&1:r^=a
            b>>=1;a<<=1
            if a&self.q:a^=self.modulus
        return r
    def power(self,a,n):
        if n<0:return self.power(self.inv(a),-n)
        r=1
        while n:
            if n&1:r=self.mul(r,a)
            a=self.mul(a,a);n>>=1
        return r
    def inv(self,a):
        if not a:raise ZeroDivisionError
        return self.power(a,self.q-2)
    def div(self,a,b):return self.mul(a,self.inv(b))


def ev(s,e,F):
    def run(n):
        if isinstance(n,ast.Name):return e[n.id]
        if isinstance(n,ast.Constant) and type(n.value)==int:return n.value&1
        if isinstance(n,ast.UnaryOp) and isinstance(n.op,ast.USub):return run(n.operand)
        if isinstance(n,ast.BinOp):
            a=run(n.left)
            if isinstance(n.op,ast.Pow):
                if not isinstance(n.right,ast.Constant):raise ValueError('Exponent')
                return F.power(a,n.right.value)
            b=run(n.right)
            if isinstance(n.op,(ast.Add,ast.Sub)):return a^b
            if isinstance(n.op,ast.Mult):return F.mul(a,b)
            if isinstance(n.op,ast.Div):return F.div(a,b)
        raise ValueError(ast.dump(n))
    return run(expression_tree(s))


@lru_cache(None)
def assignments(path,tag):
    return [(l[len(tag)+1:].split('=',1)[0].strip(),l.split('=',1)[1].lstrip('= ').strip())
            for l in (ROOT/path).read_text().splitlines() if l.startswith(tag+' ')]


def names(s):return set(re.findall('[A-Za-z_][A-Za-z_0-9]*',s))


def solve(constraints,e,F):
    e=dict(e);todo=list(constraints)
    for _ in range(len(todo)+1):
        rest=[];progress=False
        for l,r in todo:
            unknown=(names(l)|names(r))-e.keys()
            if len(unknown)==1:
                k=next(iter(unknown));roots=[]
                for v in range(F.q):
                    try:
                        if ev(l,dict(e,**{k:v}),F)==ev(r,dict(e,**{k:v}),F):roots.append(v)
                    except ZeroDivisionError:pass
                if not roots:return None
                e[k]=next((v for v in roots if v),0);progress=True
            elif not unknown:
                try:
                    if ev(l,e,F)!=ev(r,e,F):return None
                except ZeroDivisionError:return None
            else:rest.append((l,r))
        todo=rest
        if not todo:return e
        if not progress:return None
    return None


@lru_cache(None)
def coefficients(model,parameters,modulus):
    F=GF2(modulus);e=dict(parameters);result={}
    for l in (ROOT/'binary'/model/'coordinates').read_text().splitlines():
        if re.match('a[013246] =',l):
            k,v=l.split('=',1);result[k.strip()]=ev(v.strip(),e,F)
    return result


def add(P,Q,c,F):
    if P is None:return Q
    if Q is None:return P
    x,y=P;u,v=Q
    a0,a1,a2,a3,a4=(c[k] for k in ('a0','a1','a2','a3','a4'))
    if x==u:
        if v==y^F.div(F.mul(a1,x)^a3,a0):return None
        slope=F.div(F.mul(x,x)^a4^F.mul(a1,y),F.mul(a1,x)^a3)
    else:slope=F.div(y^v,x^u)
    intercept=y^F.mul(slope,x)
    X=F.mul(a0,F.mul(slope,slope))^F.mul(a1,slope)^a2^x^u
    Y=F.mul(slope,X)^intercept^F.div(F.mul(a1,X)^a3,a0)
    return X,Y


def to_weier(model,P,parameters,F):
    if model=='edwards' and P==dict(x=0,y=0):return None
    out=[ev(v,dict(parameters,**P),F) for _,v in assignments('binary/'+model+'/coordinates','toweierstrass')]
    x,y=out;c=coefficients(model,tuple(sorted(parameters.items())),F.modulus)
    check=F.mul(c['a0'],F.mul(y,y))^F.mul(c['a1'],F.mul(x,y))^F.mul(c['a3'],y)^F.power(x,3)^F.mul(c['a2'],F.mul(x,x))^F.mul(c['a4'],x)^c['a6']
    if check:raise AssertionError(('map off curve',model,P,parameters,F.q))
    return x,y


def from_weier(model,P,parameters,F):
    if P is None:
        if model=='edwards':return dict(x=0,y=0)
        raise ZeroDivisionError('infinity')
    keys=[k for k,_ in assignments('binary/'+model+'/coordinates','toweierstrass')]
    e=dict(parameters,**dict(zip(keys,P)))
    return {k:ev(v,e,F) for k,v in assignments('binary/'+model+'/coordinates','fromweierstrass')}


@lru_cache(None)
def pool(model,parameters,modulus):
    F=GF2(modulus);e=dict(parameters);relations=assignments('binary/'+model+'/coordinates','satisfying');out=[]
    for x in range(F.q):
        for y in range(F.q):
            P=dict(x=x,y=y)
            if all(ev(l,dict(e,**P),F)==ev(r,dict(e,**P),F) for l,r in relations):
                try:
                    W=to_weier(model,P,e,F)
                    if W is not None:out.append((P,W))
                except ZeroDivisionError:pass
    return out


def lift(model,coord,P,index,F,unit=False):
    x,y=P['x'],P['y'];z=1 if unit else index+2
    if model=='edwards' and coord in ('w','w-1'):return {'w'+str(index):x^y}
    if model=='edwards' and coord.startswith('wz'):e=dict(W=F.mul(x^y,z),Z=z)
    elif coord in ('affine','xy','xy-1'):e=dict(X=x,Y=y)
    elif coord=='jacobian':e=dict(X=F.mul(x,F.power(z,2)),Y=F.mul(y,F.power(z,3)),Z=z)
    elif coord=='lambda':e=dict(X=F.mul(x,z),L=F.mul(z,x^F.div(y,x)),Z=z)
    elif model=='shortw' and ('lopezdahab' in coord or coord.startswith('extended')):e=dict(X=F.mul(x,z),Y=F.mul(y,F.power(z,2)),Z=z)
    elif coord=='xz':e=dict(X=F.mul(x,z),Z=z)
    else:e=dict(X=F.mul(x,z),Y=F.mul(y,z),Z=z)
    if 'Z' in e:
        e.update(ZZ=F.mul(e['Z'],e['Z']),XZ=F.mul(e.get('X',0),e['Z']))
    return {k+str(index):v for k,v in e.items()}


def decode(model,coord,out,index,params,F):
    j=str(index)
    if model=='edwards' and coord in ('w','w-1'):return out['w'+j]
    if model=='edwards' and coord.startswith('wz'):return F.div(out['W'+j],out['Z'+j])
    if coord in ('affine','xy','xy-1'):P=dict(x=out['X'+j],y=out['Y'+j])
    else:
        z=out['Z'+j]
        if not z:raise ZeroDivisionError('affine output chart')
        x=F.div(out['X'+j],F.power(z,2) if coord=='jacobian' else z)
        if coord=='xz':return x
        if coord=='lambda':y=F.mul(x,F.div(out['L'+j]^out['X'+j],z))
        else:
            power=3 if coord=='jacobian' else 2 if model=='shortw' and ('lopezdahab' in coord or coord.startswith('extended')) else 1
            y=F.div(out['Y'+j],F.power(z,power))
        if 'ZZ'+j in out and out['ZZ'+j]!=F.mul(z,z):raise AssertionError('cached ZZ')
        if 'XZ'+j in out and out['XZ'+j]!=F.mul(out['X'+j],z):raise AssertionError('cached XZ')
        P=dict(x=x,y=y)
    return to_weier(model,P,params,F)


def execute(program,e,F):
    e=dict(e)
    for k,v in program:e[k]=ev(v,e,F)
    return e


def verify(path):
    _,model,coord,operation,_=path.split('/')
    source='binary/'+model+'/coordinates'
    params=[l.split()[1] for l in (ROOT/source).read_text().splitlines() if l.startswith('parameter ')]
    ca=assignments('binary/'+model+'/'+coord+'/variables','assume');fa=assignments(path,'assume')
    programs=[assignments(path,'compute'),[(l.split('=',1)[0].strip(),l.split('=',1)[1].strip()) for l in (ROOT/(path+'.op3')).read_text().splitlines()]]
    units={int(l[-1]) for l,r in fa if re.fullmatch('Z[123]',l) and r=='1'}
    counts=0;fields=set();excluded=0;seen_configs=set();seen_cases=set()
    for modulus in (19,37,67):
        F=GF2(modulus)
        for av in (2,3,4,5):
            e={k:(av if k in ('a2','d1','d') else 1 if k=='a6' else av^1) for k in params}
            for l,r in ca:
                if l in e and names(r)<=e.keys():e[l]=ev(r,e,F)
            e=solve(ca+[(l,r) for l,r in fa if not re.search(r'\bZ[123]\b',l+r)],e,F)
            if e is None:continue
            config_key=(modulus,tuple(sorted(e.items())))
            if config_key in seen_configs:continue
            seen_configs.add(config_key)
            base={k:e[k] for k in params}
            if model=='edwards' and (not base['d1'] or base['d2']==F.mul(base['d1'],base['d1'])^base['d1']):continue
            if model=='hessian' and F.power(base['d'],3)==1:continue
            pts=pool(model,tuple(sorted(base.items())),modulus)
            if len(pts)<6:continue
            c=coefficients(model,tuple(sorted(base.items())),modulus)
            for i in range(20):
                P,WP=pts[(7*i+1)%len(pts)];Q,WQ=pts[(11*i+3)%len(pts)]
                try:
                    env=dict(e)
                    if operation in ('diffadd','ladder'):
                        negP=(WP[0],WP[1]^F.div(F.mul(c['a1'],WP[0])^c['a3'],c['a0']))
                        D=add(WQ,negP,c,F);affD=from_weier(model,D,base,F)
                        env.update(lift(model,coord,affD,1,F,1 in units));env.update(lift(model,coord,P,2,F,2 in units));env.update(lift(model,coord,Q,3,F,3 in units))
                        expected=[(4,add(WP,WP,c,F)),(5,add(WP,WQ,c,F))] if operation=='ladder' else [(4,add(WP,WQ,c,F))]
                    else:
                        env.update(lift(model,coord,P,1,F,1 in units))
                        if operation=='addition':env.update(lift(model,coord,Q,2,F,2 in units))
                        expected=[(3,add(WP,WQ,c,F) if operation=='addition' else add(add(WP,WP,c,F),WP,c,F) if operation=='tripling' else add(WP,WP,c,F) if operation=='doubling' else WP)]
                    if solve(fa,env,F) is None:continue
                    case_key=(modulus,tuple(sorted(env.items())))
                    if case_key in seen_cases:continue
                    seen_cases.add(case_key)
                    results=[execute(pr,env,F) for pr in programs]
                    for index,target in expected:
                        if coord=='xz':wanted=from_weier(model,target,base,F)['x']
                        elif model=='edwards' and coord in ('w','w-1','wz','wz-1'):
                            targetP=from_weier(model,target,base,F);wanted=targetP['x']^targetP['y']
                        else:wanted=target
                        for out in results:
                            actual=decode(model,coord,out,index,base,F)
                            if actual!=wanted:raise AssertionError(('group mismatch',F.q,base,P,Q,index,actual,wanted))
                except (ZeroDivisionError,ValueError):excluded+=1;continue
                counts+=1;fields.add(F.q)
            if counts>=24 and len(fields)>=2:return dict(status='finite_field_checked',cases=counts,fields=sorted(fields),chart_or_denominator_exclusions=excluded)
    return dict(status='insufficient_coverage',cases=counts,fields=sorted(fields),chart_or_denominator_exclusions=excluded)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'research/binary-verification.json');args=p.parse_args()
    records=[]
    for r in json.loads((ROOT/'research/binary-import-manifest.json').read_text())['files']:
        if not r['operation_code']:continue
        try:result=verify(r['path'])
        except Exception as error:result=dict(status='failed',error=repr(error))
        if result['status']!='finite_field_checked':print(r['path'],result,flush=True)
        records.append(dict(path=r['path'],metadata_sha256=hashlib.sha256((ROOT/r['path']).read_bytes()).hexdigest(),op3_sha256=hashlib.sha256((ROOT/(r['path']+'.op3')).read_bytes()).hexdigest(),**result))
    args.output.write_text(json.dumps(dict(method='Exact polynomial-basis GF(2^m) arithmetic over multiple extension degrees, checked against general Weierstrass group law. Affine chart coverage only; no completeness or symbolic proof claimed.',formulas=records),indent=2)+'\n')
    counts=Counter(r['status'] for r in records);print(dict(counts),flush=True)
    if set(counts)!={'finite_field_checked'}:raise SystemExit(1)


if __name__=='__main__':main()
