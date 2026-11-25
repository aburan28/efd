#!/usr/bin/env python3
"""Inspect the atlas without evaluating formula expressions or adding dependencies."""
import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--curve', help='Curve ID or alias (e.g. p256, curve25519)')
    parser.add_argument('--json', action='store_true', help='Print the complete catalog, or selected curve')
    args = parser.parse_args()
    catalog = json.loads((ROOT / 'transformations/catalog.json').read_text())
    curve = None
    if args.curve:
        curve = next((c for c in catalog['curves'] if args.curve.lower() in
                      [c['id'].lower(), *[a.lower() for a in c.get('aliases', [])]]), None)
        if curve is None:
            parser.error('Unknown curve; choose ' + ', '.join(c['id'] for c in catalog['curves']))
    if args.json:
        print(json.dumps(curve or catalog, indent=2))
        return
    if curve:
        print(curve['id'] + ': ' + ', '.join(curve['base_field_models']))
        print('Coordinate candidates: ' + ', '.join(curve['recommended_coordinates']))
        for condition in curve['conditions']:
            print('  ' + condition)
        return
    for edge in catalog['transformations']:
        print(f"{edge['id']}: {edge['source']} <-> {edge['target']}")
        print('  Requires: ' + '; '.join(edge['requires']))
    print('\nCoordinate representations (see native variables files for equations):')
    for coord in catalog['coordinates']:
        print(f"  {coord['id']}: {coord['kind']}; {coord['affine_chart']}")
    print('\nAdditional native formulas:')
    for formula in catalog.get('additional_formula_files', []):
        print('  ' + formula['path'])
    print('\nDegree-two isogeny and dual: transformations/two-isogeny.json')


if __name__ == '__main__':
    main()
