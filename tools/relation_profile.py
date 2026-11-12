#!/usr/bin/env python3
"""Measure a three-point summation relation; this is not an ECDLP solver."""
import argparse
import json
from pathlib import Path
import sympy as S

ROOT = Path(__file__).resolve().parents[1]


def summation_polynomial(x1, x2, x3, a, b):
    return ((x1-x2)**2*x3**2 - 2*((x1+x2)*(x1*x2+a)+2*b)*x3
            + (x1*x2-a)**2 - 4*b*(x1+x2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--curve', choices=['curve25519', 'p256', 'secp256k1'], default='curve25519')
    args = parser.parse_args()
    cat = json.loads((ROOT / 'transformations/catalog.json').read_text())
    curve = next(c for c in cat['curves'] if c['id'] == args.curve)
    p = int(curve['p'])
    a, b = (int(curve['parameters'][k]) % p for k in ('a', 'b'))
    xs = S.symbols('x1:4')
    x1,x2,x3 = xs
    if curve['model'] == 'montgom':
        A,B = a,b
        a,b = ((1-A*A*pow(3,-1,p))*pow(B,-2,p) % p,
               (2*A**3*pow(27,-1,p)-A*pow(3,-1,p))*pow(B,-3,p) % p)
    f3 = summation_polynomial(x1, x2, x3, a, b)
    variants = [('shortw_x', f3)]
    if curve['model'] == 'montgom':
        shift = A*pow(3,-1,p) % p
        variants.append(('montgomery_u', f3.subs({x:(x+shift)*pow(B,-1,p) for x in xs}, simultaneous=True)))
        variants.append(('edwards_y', f3.subs({x:((1+x)/(1-x)+shift)*pow(B,-1,p) for x in xs}, simultaneous=True)))
    results=[]
    for name, relation in variants:
        numerator = S.cancel(relation).as_numer_denom()[0]
        poly = S.Poly(numerator, *xs, modulus=p)
        results.append(dict(coordinates=name, total_degree=poly.total_degree(),
                            multidegree=[poly.degree(x) for x in xs], terms=len(poly.terms())))
    print(json.dumps(dict(curve=curve['id'], relation='Semaev f3; sign choices for full points are implicit',
                         measurements=results,
                         limitation='Equation metrics only; no factor base, relation probability, solver timing or DLP speedup measured. Edwards chart excludes y_i=1.'), indent=2))


if __name__ == '__main__':
    main()
