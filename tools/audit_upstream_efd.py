#!/usr/bin/env python3
"""Inventory upstream EFD links without importing or executing upstream formulas."""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from datetime import date
from html.parser import HTMLParser
import json
from pathlib import Path
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[1]
BASE = 'https://www.hyperelliptic.org/EFD/g1p/'


class Links(HTMLParser):
    def __init__(self, html):
        super().__init__()
        self.links = set()
        self.feed(html)

    def handle_starttag(self, tag, attrs):
        href = dict(attrs).get('href', '')
        if tag == 'a': self.links.add(href)


def fetch(name):
    with urlopen(BASE+name, timeout=30) as response:
        return name, Links(response.read().decode('utf-8'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    _, index = fetch('index.html')
    pages = sorted(h for h in index.links if h.startswith('auto-') and h.endswith('.html'))
    references = set()
    # Four workers bound the load on the upstream server; every failed page aborts the audit.
    with ThreadPoolExecutor(max_workers=4) as pool:
        for _, parsed in pool.map(fetch, pages):
            references.update(h[5:] for h in parsed.links if h.startswith('data/'))
    formulas = sorted(h for h in references if len(h.split('/')) == 4)
    models = sorted(h for h in references if h.endswith('/coordinates'))
    coordinates = sorted(h for h in references if h.endswith('/variables'))
    missing = [h for h in formulas if not (ROOT/h).is_file()]
    native = sorted(str(p.relative_to(ROOT)) for model in ('shortw','montgom','twisted','edwards')
                    for p in (ROOT/model).glob('*/*/*') if p.is_file() and p.suffix != '.op3')
    result = dict(checked_on=date.today().isoformat(), source=BASE+'index.html',
                  comparison='Exact relative paths; presence does not prove algebraic equality or source-version equality. Missing entries may include formulas deliberately removed as broken by the export.',
                  status='Discovery inventory only; no upstream formulas imported or verified.',
                  url_templates=dict(metadata=BASE+'data/{path}',operation_code=BASE+'auto-code/{path}.op3'),
                  counts=dict(model_families=len(models),coordinate_representations=len(coordinates),
                              formula_links=len(formulas),local_formula_files=len(native),
                              missing_formula_paths=len(missing)),
                  missing_models=[h for h in models if not (ROOT/h).is_file()],
                  missing_coordinates=[h for h in coordinates if not (ROOT/h).is_file()],
                  missing_formulas_by_model=dict(sorted(Counter(h.split('/')[0] for h in missing).items())),
                  missing_formulas=missing,
                  local_formulas_not_linked_upstream=sorted(set(native)-set(formulas)))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('counts','missing_models','missing_formulas_by_model')},indent=2))


if __name__ == '__main__': main()
