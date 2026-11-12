"""Check new native formula code against independent affine group arithmetic."""
import ast
import json
from pathlib import Path
import unittest
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
CAT = json.loads((ROOT/'transformations/catalog.json').read_text())


def add(P, Q, p, a2, a4):
    if P is None: return Q
    if Q is None: return P
    x,y = P; u,v = Q
    if x == u and (y+v) % p == 0: return None
    slope = ((3*x*x+2*a2*x+a4)*pow(2*y,-1,p) if P==Q
             else (v-y)*pow(u-x,-1,p)) % p
    X = (slope*slope-a2-x-u) % p
    return X, (slope*(x-X)-y) % p


def evaluate(path, env, p):
    env = dict(env)
    for line in path.read_text().splitlines():
        if path.suffix != '.op3':
            if not line.startswith('compute '): continue
            line = line[len('compute '):]
        lhs, rhs = line.split('=',1)
        # Only checked-in, trusted arithmetic expressions, never user input.
        env[lhs.strip()] = int(S.sympify(rhs.replace('^','**'), locals=env)) % p
    return env


def same_x(test, result, X, Z, P, p):
    test.assertNotEqual((result[X],result[Z]), (0,0), 'undefined projective output')
    if P is None: test.assertEqual(result[Z], 0)
    else: test.assertEqual((result[X]-P[0]*result[Z]) % p, 0)


class AdditionalFormulas(unittest.TestCase):
    def test_op3_structure_and_metadata_agree(self):
        for rec in CAT['additional_formula_files']:
            path = ROOT/rec['path']
            for line in path.with_suffix('.op3').read_text().splitlines():
                tree = ast.parse(line.replace('^','**')).body[0]
                self.assertIsInstance(tree, ast.Assign)
                self.assertIsInstance(tree.value, ast.BinOp)
                self.assertIsInstance(tree.value.left, (ast.Name, ast.Constant))
                self.assertIsInstance(tree.value.right, (ast.Name, ast.Constant))
            compute = [l[8:] for l in path.read_text().splitlines() if l.startswith('compute ')]
            self.assertEqual(compute, path.with_suffix('.op3').read_text().splitlines())

    def test_montgomery_doubling_and_ladders(self):
        # Generic projective rescalings and several nonsingular toy curves.
        # Absorb b using x_SW=b*x_M, y_SW=b^2*y_M, giving a4=b^2 and a2=a*b.
        # This covers b!=1 without a field square root.
        for p,a,b in [(101,5,1),(103,7,3),(107,9,2)]:
            points=[(x,y) for x in range(1,p) for y in range(1,p)
                    if (b*y*y-x**3-a*x*x-x) % p == 0][:8]
            def group(P,Q):
                conv=lambda P: None if P is None else (b*P[0]%p,b*b*P[1]%p)
                R=add(conv(P),conv(Q),p,a*b%p,b*b%p)
                return None if R is None else (R[0]*pow(b,-1,p)%p,R[1]*pow(b*b,-1,p)%p)
            for P in points:
                for Q in points:
                    if P[0] == Q[0]: continue  # Difference O or 2-torsion needs a separate ladder case.
                    diff=group(Q,(P[0],-P[1]%p))
                    if diff is None or diff[0] == 0: continue
                    for rec in CAT['additional_formula_files']:
                        if not rec['path'].startswith('montgom/'): continue
                        path=ROOT/rec['path']
                        mixed=path.name.startswith('mladd')
                        z1,z2,z3=(1 if mixed else 7),11,13
                        env=dict(a=a,a24=(a-2)*pow(4,-1,p)%p,
                                 X1=diff[0]*z1%p,Z1=z1,X2=P[0]*z2%p,Z2=z2,
                                 X3=Q[0]*z3%p,Z3=z3)
                        if rec['kind']=='doubling': env.update(X1=P[0]*z2%p,Z1=z2)
                        for formula in (path,path.with_suffix('.op3')):
                            result=evaluate(formula,env,p)
                            with self.subTest(formula=str(formula),p=p):
                                if rec['kind']=='doubling': same_x(self,result,'X3','Z3',group(P,P),p)
                                else:
                                    same_x(self,result,'X4','Z4',group(P,P),p)
                                    same_x(self,result,'X5','Z5',group(P,Q),p)

    def test_shortw_tripling_all_toy_points_and_infinity(self):
        path=ROOT/'shortw/xz/tripling/tpl-division'
        for p,a,b in [(101,2,3),(103,-3,5),(107,0,7)]:
            points=[None]+[(x,y) for x in range(p) for y in range(p)
                          if (y*y-x**3-a*x-b) % p == 0]
            # y and -y have the same quotient coordinate and output.
            seen=set()
            for P in points:
                key=None if P is None else P[0]
                if key in seen: continue
                seen.add(key)
                for scale in (1,7):
                    X,Z=(scale,0) if P is None else (P[0]*scale%p,scale)
                    env=dict(a=a,b=b,X1=X,Z1=Z)
                    expected=add(add(P,P,p,0,a),P,p,0,a)
                    for formula in (path,path.with_suffix('.op3')):
                        result=evaluate(formula,env,p)
                        with self.subTest(formula=str(formula),p=p,point=P):
                            same_x(self,result,'X3','Z3',expected,p)

    def test_two_isogeny_symbolic_identity_and_dual(self):
        iso=json.loads((ROOT/'transformations/two-isogeny.json').read_text())
        x,y,A,B=S.symbols('x y A B'); names=dict(x=x,y=y,A=A,B=B)
        f={k:S.sympify(v,locals=names) for k,v in iso['forward'].items()}
        dual={k:S.sympify(iso['dual'][k],locals=names) for k in ('x','y')}
        relation=y*y-x**3-A*x*x-B*x
        target=f['y']**2-f['x']**3+2*A*f['x']**2-(A*A-4*B)*f['x']
        rem=S.Poly(S.cancel(target).as_numer_denom()[0],y).rem(S.Poly(relation,y)).as_expr()
        self.assertEqual(S.cancel(rem),0)
        slope=(3*x*x+2*A*x+B)/(2*y)
        double_x=slope*slope-A-2*x
        double_y=slope*(x-double_x)-y
        for key,expected in [('x',double_x),('y',double_y)]:
            result=dual[key].subs({x:f['x'],y:f['y']},simultaneous=True)-expected
            num=S.cancel(result).as_numer_denom()[0]
            rem=S.Poly(num,y).rem(S.Poly(relation,y)).as_expr()
            self.assertEqual(S.cancel(rem),0)

    def test_two_isogeny_homomorphism_with_kernel(self):
        for p,A,B in [(101,5,1),(103,7,2)]:
            points=[None]+[(x,y) for x in range(p) for y in range(p)
                          if (y*y-x**3-A*x*x-B*x) % p == 0][:16]
            def phi(P):
                if P is None or P[0]==0: return None
                x,y=P
                return ((x+A+B*pow(x,-1,p))%p,y*(1-B*pow(x,-2,p))%p)
            for P in points:
                for Q in points:
                    self.assertEqual(phi(add(P,Q,p,A,B)),
                                     add(phi(P),phi(Q),p,-2*A,(A*A-4*B)%p))


if __name__=='__main__': unittest.main()
