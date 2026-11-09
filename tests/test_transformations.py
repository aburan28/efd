"""Exact identities for trusted local catalog expressions; no numerical floats."""
import json
from pathlib import Path
import unittest
import sympy as S
from tools.relation_profile import summation_polynomial

ROOT = Path(__file__).resolve().parents[1]
CAT = json.loads((ROOT / 'transformations/catalog.json').read_text())
NAMES = 'x y a b c d r s t X Y Z T ZZ ZZZ u v weierx weiery L'.split()
SYMS = {n: S.Symbol(n) for n in NAMES}
x, y = SYMS['x'], SYMS['y']


def expr(value):
    return S.sympify(value.replace('^', '**'), locals=SYMS)


def equations(model, params=None, point=None):
    params = params or {k: SYMS[k] for k in 'abcd'}
    px, py = point or (x, y)
    a, b, c, d = (params.get(k, SYMS[k]) for k in 'abcd')
    return {'shortw': py**2-px**3-a*px-b,
            'montgom': b*py**2-px**3-a*px**2-px,
            'twisted': a*px**2+py**2-1-d*px**2*py**2,
            'edwards': px**2+py**2-c**2*(1+d*px**2*py**2)}[model]


def assignments(path, tag):
    result = {}
    for line in path.read_text().splitlines():
        if line.startswith(tag + ' '):
            lhs, rhs = line[len(tag)+1:].split('=', 1)
            result[expr(lhs.strip())] = expr(rhs.strip())
    return result


def resolve(value, definitions):
    for _ in range(len(definitions) + 1):
        new = value.subs(definitions, simultaneous=True)
        if new == value:
            return new
        value = new
    return value


class Transformations(unittest.TestCase):
    def test_summation_relation_on_toy_curve_additions(self):
        p, a, b = 101, 2, 3
        points = [(u,v) for u in range(p) for v in range(p)
                  if (v*v-u**3-a*u-b) % p == 0]
        for u,v in points:
            for w,z in points:
                if u == w:
                    if (v+z) % p == 0: continue
                    slope = (3*u*u+a)*pow(2*v,-1,p) % p
                else:
                    slope = (z-v)*pow(w-u,-1,p) % p
                result_x = (slope*slope-u-w) % p
                self.assertEqual(summation_polynomial(u,w,result_x,a,b) % p, 0)

    def test_native_weierstrass_maps(self):
        for model in CAT['models']:
            with self.subTest(model=model):
                path = ROOT / model / 'coordinates'
                forward = assignments(path, 'toweierstrass')
                inverse = assignments(path, 'fromweierstrass')
                invariants = {}
                for line in path.read_text().splitlines():
                    if line.startswith(('a0 =', 'a1 =', 'a2 =', 'a3 =', 'a4 =', 'a6 =')):
                        key,value = line.split('=',1)
                        invariants[key.strip()] = expr(value.strip())
                px = forward.get(SYMS['u'], forward.get(SYMS['weierx']))
                py = forward.get(SYMS['v'], forward.get(SYMS['weiery']))
                target = (invariants['a0']*py**2 + invariants['a1']*px*py
                          + invariants['a3']*py - px**3 - invariants['a2']*px**2
                          - invariants['a4']*px - invariants['a6'])
                num = S.cancel(target).as_numer_denom()[0]
                rem = S.Poly(num,y).rem(S.Poly(equations(model),y)).as_expr()
                self.assertEqual(S.cancel(rem), 0)
                for variable, formula in inverse.items():
                    self.assertEqual(S.cancel(formula.subs(forward, simultaneous=True)-variable), 0)

    def test_catalog_links_and_coverage(self):
        self.assertEqual({c['path'] for c in CAT['coordinates']},
                         {str(p.relative_to(ROOT)) for p in ROOT.glob('*/*/variables')})
        self.assertEqual(len({m['id'] for m in CAT['transformations']}), len(CAT['transformations']))
        for record in [*CAT['models'].values(), *CAT['coordinates']]:
            self.assertTrue((ROOT / record['path']).is_file())
        for record in [*CAT['transformations'], *CAT['curves'], *CAT['coordinates']]:
            self.assertTrue(all(s in CAT['sources'] for s in record['sources']))
        ids = {c['id'] for c in CAT['coordinates']}
        for curve in CAT['curves']:
            self.assertTrue(set(curve['recommended_coordinates']) <= ids)

    def test_curve_equations_and_both_inverse_compositions(self):
        for m in CAT['transformations']:
            with self.subTest(map=m['id']):
                subs = {expr(k): expr(v) for k, v in m['source_substitutions'].items()}
                f = {expr(k): expr(v) for k, v in m['forward'].items()}
                inv = {expr(k): expr(v) for k, v in m['inverse'].items()}
                params = {k: expr(v) for k, v in m['target_parameters'].items()}
                relation = equations(m['source']).subs(subs, simultaneous=True)
                target = equations(m['target'], params, (f[x], f[y])).subs(subs, simultaneous=True)
                numerator = S.cancel(target).as_numer_denom()[0]
                remainder = S.Poly(numerator, y).rem(S.Poly(relation, y)).as_expr()
                self.assertEqual(S.cancel(remainder), 0)
                for outer, inner in [(inv, f), (f, inv)]:
                    for variable in (x, y):
                        result = outer[variable].subs(inner, simultaneous=True)
                        self.assertEqual(S.cancel((result-variable).subs(subs, simultaneous=True)), 0)

    def test_native_coordinate_roundtrips(self):
        for coord in CAT['coordinates']:
            path = ROOT / coord['path']
            back = assignments(path, 'toaffine')
            if not back:  # Quotient coordinates intentionally have no full inverse.
                continue
            with self.subTest(coordinates=coord['id']):
                definitions = assignments(path, 'satisfying')
                lift = assignments(path, 'tosystem')
                for variable, formula in back.items():
                    expanded = resolve(formula, {k:v for k,v in definitions.items() if k not in (x,y)})
                    self.assertEqual(S.cancel(expanded.subs(lift, simultaneous=True)-variable), 0)
                    self.assertEqual(S.cancel(expanded-resolve(definitions[variable], {k:v for k,v in definitions.items() if k not in (x,y)})), 0)

    def test_homogeneous_weights(self):
        L = SYMS['L']
        for coord in CAT['coordinates']:
            path = ROOT / coord['path']
            weights = assignments(path, 'homogweight')
            scale = {var: var*L**weight for var, weight in weights.items()}
            relations = assignments(path, 'satisfying')
            for lhs, rhs in relations.items():
                # Expand cached ZZ/ZZZ to base variables before checking degree.
                extra = {k:v for k,v in relations.items() if k not in (x,y) and k not in weights}
                value = resolve(lhs-rhs, extra).subs(scale, simultaneous=True)
                numerator = S.cancel(value).as_numer_denom()[0]
                with self.subTest(coordinates=coord['id'], relation=str(lhs)):
                    self.assertLessEqual(len(S.Poly(numerator, L).terms()), 1)

    def test_real_curves_and_scalar_orders(self):
        for curve in CAT['curves']:
            with self.subTest(curve=curve['id']):
                p, n = int(curve['p']), int(curve['subgroup_order'])
                pars = {k:int(v) % p for k,v in curve['parameters'].items()}
                P = tuple(int(curve['base_point'][k]) for k in ('x','y'))
                a, b = pars['a'], pars['b']
                if curve['model'] == 'montgom':
                    self.assertEqual((b*P[1]**2-P[0]**3-a*P[0]**2-P[0]) % p, 0)
                    # Map to short Weierstrass for a single independent group-law check.
                    P = ((P[0]+a*pow(3,-1,p))*pow(b,-1,p) % p, P[1]*pow(b,-1,p) % p)
                    a, b = ((1-a*a*pow(3,-1,p))*pow(b,-2,p) % p,
                            (2*a**3*pow(27,-1,p)-a*pow(3,-1,p))*pow(b,-3,p) % p)
                self.assertEqual((P[1]**2-P[0]**3-a*P[0]-b) % p, 0)
                self.assertNotEqual((4*a**3+27*b*b) % p, 0)
                def add(P, Q):
                    if P is None: return Q
                    if Q is None: return P
                    u,v = P; w,z = Q
                    if u == w and (v+z) % p == 0: return None
                    slope = ((3*u*u+a)*pow(2*v,-1,p) if P==Q else (z-v)*pow(w-u,-1,p)) % p
                    X = (slope*slope-u-w) % p
                    return X, (slope*(u-X)-v) % p
                R, Q, scalar = None, P, n
                while scalar:
                    if scalar & 1: R = add(R,Q)
                    Q = add(Q,Q); scalar >>= 1
                self.assertIsNone(R)
                if curve['cofactor'] == 1:
                    # Exact absence of rational 2-torsion, independent of trusting h=1.
                    modulus = S.Poly(x**3+a*x+b, x, modulus=p)
                    xp = S.Poly(1,x,modulus=p)
                    base = S.Poly(x,x,modulus=p)
                    power = p
                    while power:
                        if power & 1: xp = (xp*base).rem(modulus)
                        base = (base*base).rem(modulus); power >>= 1
                    self.assertEqual(S.gcd(modulus, xp-S.Poly(x,x,modulus=p)).degree(), 0)

    def test_rfc7748_normalization(self):
        curve = CAT['curves'][0]; p = int(curve['p'])
        d = -121665*pow(121666,-1,p) % p
        self.assertEqual(d, 37095705934669439343138083508754565189542113879843219016388785533085940283555)
        self.assertEqual(pow(-486664 % p, (p-1)//2, p), 1)
        # RFC's Edwards base point maps to the Montgomery u=9, with a consistent root sign.
        ex = 15112221349535400772501151409588531511454012693041857206046113283949847762202
        ey = 46316835694926478169428394003475163141307993866256225615783033603165251855960
        u, v = int(curve['base_point']['x']), int(curve['base_point']['y'])
        s = ex*v*pow(u,-1,p) % p
        self.assertEqual((s*s+486664) % p, 0)
        self.assertEqual((1+ey)*pow(1-ey,-1,p) % p, u)
        self.assertEqual((-ex*ex+ey*ey-1-d*ex*ex*ey*ey) % p, 0)


if __name__ == '__main__':
    unittest.main()
