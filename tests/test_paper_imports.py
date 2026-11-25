"""Independent group checks for paper-derived joint states, quotients and recovery."""
import json
from pathlib import Path
import unittest
from tools.efd_format import evaluate

ROOT=Path(__file__).resolve().parents[1]
REC=json.loads((ROOT/'research/paper-imports.json').read_text())['formulas']


def run(path,env,p):
    env=dict(env)
    for line in (ROOT/path).read_text().splitlines():
        if path.endswith('.op3'):assignment=line
        elif line.startswith('compute '):assignment=line[8:]
        else:continue
        lhs,rhs=assignment.split('=',1);env[lhs.strip()]=evaluate(rhs.strip(),env,p)
    return env


def sw_add(P,Q,p,a):
    if P is None:return Q
    if Q is None:return P
    x,y=P;u,v=Q
    if x==u and (y+v)%p==0:return None
    slope=((3*x*x+a)*pow(2*y,-1,p) if P==Q else (v-y)*pow(u-x,-1,p))%p
    X=(slope*slope-x-u)%p
    return X,(slope*(x-X)-y)%p


def mont_add(P,Q,p,a,b):
    if P is None:return Q
    if Q is None:return P
    x,y=P;u,v=Q
    if x==u and (y+v)%p==0:return None
    slope=((3*x*x+2*a*x+1)*pow(2*b*y,-1,p) if P==Q else (v-y)*pow(u-x,-1,p))%p
    X=(b*slope*slope-a-x-u)%p
    return X,(slope*(x-X)-y)%p


def ed_add(P,Q,p,a,d):
    x,y=P;u,v=Q;t=d*x*y*u*v%p
    return ((x*v+y*u)*pow(1+t,-1,p)%p,(y*v-a*x*u)*pow(1-t,-1,p)%p)


def mont_lift(P,p):
    if P is None:return (0,1,0,0)
    x,y=P
    return (0,0,1,0) if x==0 else (x,y,y*pow(x,-1,p)%p,1)


class PaperImports(unittest.TestCase):
    def programs(self,prefix):
        return [r['path']+suffix for r in REC if r['path'].startswith(prefix) for suffix in ('','.op3')]

    def test_extended_montgomery_complete_conditions_and_linear_embedding(self):
        for p in (101,109,157):
            a=next(a for a in range(3,p) if pow(a+2,(p-1)//2,p)==1 and pow(a-2,(p-1)//2,p)==p-1)
            points=[None]+[(x,y) for x in range(p) for y in range(p) if (y*y-x**3-a*x*x-x)%p==0]
            for i,P in enumerate(points):
                Q=points[(i*7+3)%len(points)]
                for scale1,scale2 in [(1,1),(3,5)]:
                    env=dict(a=a,b=1,ta=(a+2)%p,td=(a-2)%p)
                    for index,point,scale in [(1,P,scale1),(2,Q,scale2)]:
                        env.update({k+str(index):v*scale%p for k,v in zip(('U','V','S','W'),mont_lift(point,p))})
                    expected=mont_lift(mont_add(P,Q,p,a,1),p)
                    for path in self.programs('montgom/extended/'):
                        result=run(path,env,p);actual=tuple(result[k+'3'] for k in ('U','V','S','W'))
                        self.assertNotEqual(actual,(0,0,0,0))
                        for j in range(4):
                            for k in range(4):self.assertEqual((actual[j]*expected[k]-actual[k]*expected[j])%p,0)
                        U,V,S,W=actual
                        self.assertEqual((S*U-V*W)%p,0)
                        self.assertEqual((S*V-U*U-a*U*W-W*W)%p,0)

    def test_joint_coz_addition_and_updated_denominators(self):
        for p,a,b in [(101,2,3),(109,-3,7),(157,0,7)]:
            points=[(x,y) for x in range(p) for y in range(p) if (y*y-x**3-a*x-b)%p==0]
            for i,P in enumerate(points[:35]):
                Q=points[(i*7+11)%len(points)]
                if P[0]==Q[0]:continue
                z=7;env=dict(a=a,b=b,Z1=z,Z2=z,X1=P[0]*z*z%p,Y1=P[1]*z**3%p,X2=Q[0]*z*z%p,Y2=Q[1]*z**3%p)
                for prefix in ('shortw/jacobian/conjugate/','shortw/jacobian/updateadd/'):
                    expected3=sw_add(P,Q,p,a)
                    expected4=sw_add(P,(Q[0],-Q[1]%p),p,a) if 'conjugate' in prefix else P
                    for path in self.programs(prefix):
                        out=run(path,env,p)
                        self.assertEqual(out['Z3'],out['Z4'])
                        for index,expected in [(3,expected3),(4,expected4)]:
                            Z=out['Z'+str(index)];self.assertNotEqual(Z,0)
                            self.assertEqual(out['X'+str(index)]*pow(Z,-2,p)%p,expected[0])
                            self.assertEqual(out['Y'+str(index)]*pow(Z,-3,p)%p,expected[1])

    def test_recovery_and_weierstrass_relation(self):
        for model in ('shortw','montgom'):
            for p,a,b in [(101,5,1),(109,7,3)]:
                points=[(x,y) for x in range(p) for y in range(1,p)
                        if ((y*y-x**3-a*x-b) if model=='shortw' else (b*y*y-x**3-a*x*x-x))%p==0]
                for i,P in enumerate(points[:25]):
                    Q=points[(i*11+7)%len(points)]
                    R=sw_add(P,Q,p,a) if model=='shortw' else mont_add(P,Q,p,a,b)
                    if R is None:continue
                    env=dict(a=a,b=b,X1=P[0],Y1=P[1],Z1=1,X2=Q[0]*5%p,Z2=5,X3=R[0]*7%p,Z3=7)
                    for path in self.programs(model+'/projective/recovery/'):
                        out=run(path,env,p);z=out['Z4'];self.assertNotEqual(z,0)
                        self.assertEqual(out['X4']*pow(z,-1,p)%p,Q[0])
                        self.assertEqual(out['Y4']*pow(z,-1,p)%p,Q[1])
                    for path in self.programs(model+'/xz/relation/'):
                        self.assertEqual(run(path,env,p)['R'],0)
                        self.assertTrue(any(run(path,dict(env,X3=t),p)['R']!=0 for t in range(5)))

    def test_edwards_invariants_and_summation_relations(self):
        for p in (101,109,157):
            a=p-1;d=next(d for d in range(2,p) if pow(d,(p-1)//2,p)==p-1)
            points=[(x,y) for x in range(p) for y in range(p) if (a*x*x+y*y-1-d*x*x*y*y)%p==0]
            def w(P):return d*P[0]**2*P[1]**2%p
            checked=0
            for i,P in enumerate(points):
                Q=points[(i*7+3)%len(points)];R=ed_add(P,Q,p,a,d);D=ed_add(Q,(-P[0]%p,P[1]),p,a,d)
                env=dict(a=a,d=d,Y1=P[1]*3%p,Z1=3,Y2=Q[1]*5%p,Z2=5,Y3=R[1]*7%p,Z3=7)
                for path in self.programs('twisted/yz/relation/'):self.assertEqual(run(path,env,p)['R'],0)
                if w(P)**2%p!=1:
                    env=dict(a=a,d=d,e=4*a*pow(d,-1,p)%p,W1=w(P)*3%p,Z1=3)
                    expected=w(ed_add(P,P,p,a,d))
                    for path in self.programs('twisted/w/doubling/'):
                        out=run(path,env,p);self.assertNotEqual((out['W3'],out['Z3']),(0,0));self.assertEqual((out['W3']-expected*out['Z3'])%p,0)
                if w(D)==0 or w(P)*w(Q)%p==1:continue
                env=dict(a=a,d=d,W1=w(D)*3%p,Z1=3,W2=w(P)*5%p,Z2=5,W3=w(Q)*7%p,Z3=7)
                for path in self.programs('twisted/w/diffadd/'):
                    out=run(path,env,p);self.assertNotEqual((out['W4'],out['Z4']),(0,0));self.assertEqual((out['W4']-w(R)*out['Z4'])%p,0)
                checked+=1
            self.assertGreater(checked,10)

    def test_linear_maps_symbolically(self):
        import sympy as S
        X,Y,T,Z,a,d=S.symbols('X Y T Z a d')
        U,V,C,W=X+T,Z+Y,Z-Y,X-T
        self.assertEqual(S.expand(C*U-V*W),2*(Z*T-X*Y))
        self.assertEqual(S.cancel((a-d)*(4/(a-d)*C*V-U**2-2*(a+d)/(a-d)*U*W-W**2)
                                  +4*(a*X**2+Y**2-Z**2-d*T**2)),0)
        self.assertEqual(tuple(S.expand(v) for v in (U+W,V-C,U-W,V+C)),(2*X,2*Y,2*T,2*Z))

    def test_hamburg_coupled_state_against_affine_triples(self):
        checked=0
        for p,a,b in [(101,-3,5),(109,2,3),(157,7,11)]:
            points=[(x,y) for x in range(p) for y in range(p) if (y*y-x**3-a*x-b)%p==0]
            for i,P in enumerate(points[:35]):
                Q=points[(i*11+7)%len(points)]
                R=sw_add(P,Q,p,a)
                if R is None or Q[0]==R[0]:continue
                S=sw_add(Q,R,p,a);T=sw_add(R,R,p,a)
                if S is None or T is None or S[0]==T[0]:continue
                m=(Q[1]+R[1])*pow(Q[0]-R[0],-1,p)%p
                for z in (1,7):
                    env=dict(XQP=(Q[0]-P[0])*z*z%p,XRP=(R[0]-P[0])*z*z%p,M=m*z%p,YP=2*P[1]*z**3%p,Z=z)
                    for path in self.programs('shortw/jacobian/coupledstep/'):
                        out=run(path,env,p);zz=out['Zout']
                        if zz==0:continue
                        self.assertEqual(out['XSP'],(S[0]-P[0])*zz*zz%p)
                        self.assertEqual(out['XTP'],(T[0]-P[0])*zz*zz%p)
                        self.assertEqual(out['YPout'],2*P[1]*zz**3%p)
                        self.assertEqual(out['Mout'],(S[1]+T[1])*pow(S[0]-T[0],-1,p)*zz%p)
                        checked+=1
        self.assertGreater(checked,100)

    def test_general_weierstrass_chord(self):
        for p in (101,109):
            pars=dict(a1=2,a2=3,a3=5,a4=7,a6=11)
            points=[(x,y) for x in range(p) for y in range(p)
                    if (y*y+2*x*y+5*y-x**3-3*x*x-7*x-11)%p==0]
            for i,P in enumerate(points[:30]):
                Q=points[(i*7+3)%len(points)]
                if P[0]==Q[0]:continue
                slope=(Q[1]-P[1])*pow(Q[0]-P[0],-1,p)%p
                intercept=(P[1]*Q[0]-Q[1]*P[0])*pow(Q[0]-P[0],-1,p)%p
                x=(slope*slope+2*slope-3-P[0]-Q[0])%p
                expected=(x,(-(slope+2)*x-intercept-5)%p)
                env=dict(pars,X1=P[0]*3%p,Y1=P[1]*3%p,Z1=3,X2=Q[0]*7%p,Y2=Q[1]*7%p,Z2=7)
                for path in self.programs('generalw/projective/addition/'):
                    out=run(path,env,p);z=out['Z3'];self.assertNotEqual(z,0)
                    self.assertEqual((out['X3']*pow(z,-1,p)%p,out['Y3']*pow(z,-1,p)%p),expected)


if __name__=='__main__':unittest.main()
