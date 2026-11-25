# Coverage and verification

Snapshot: 2026-11-12. The [unified index](../database/index.json) contains 639
native formula entries, 16 namespaced equation families and 59 coordinate
representations. The [search tool](../tools/search_formulas.py) filters by field,
model, operation, source and verification status. These are entry counts;
rescalings, specializations and algebraically equivalent formulas are retained
as useful implementations, rather than counted as independent discoveries.

## What is included

| Section | Models | Coordinate representations | Formula entries |
| --- | ---: | ---: | ---: |
| Prime/large-characteristic | 13 | 40 | 504 |
| Ordinary binary | 3 | 19 | 135 |
| Total | 16 | 59 | 639 |

Prime-field equations include short and general Weierstrass, Montgomery,
Edwards/twisted Edwards, both Doche–Icart–Kohel families, Hessian/twisted Hessian,
Jacobi quartic/intersection, Legendre and generalized Huff. Legendre and
Huff contribute equation/map/chart definitions; there are no optimized native
arithmetic entries for them yet. Their inclusion must not be mistaken for new
prime-field models available on P-256: rational torsion and square conditions
still apply.

All 471 linked large-characteristic EFD formula paths are present. The
[latest live audit](../research/upstream-efd-current.json) reports no missing
models, coordinate systems or formulas from that section. The separate binary
import covers all 135 formula links, 19 representations and three families
found in the ordinary-binary EFD. The original four-family export's identifiers
and operation code are preserved, apart from the coordinate/map corrections
already described in the PR.

The new imports include:

- 200 previously absent prime-field EFD entries: 194 omitted-family entries and
  six Chudnovsky additions. The old broken operation code is replaced by fresh
  arithmetic-AST compilation of the source metadata, not restored as broken text.
- 135 binary entries: binary Edwards, Hessian and ordinary Weierstrass,
  including affine, projective, Jacobian, Lopez–Dahab, extended, lambda and
  quotient representations.
- Extended Montgomery full-point and standard projective representations;
  Edwards y and `w=d*x^2*y^2` quotient representations.
- Twelve additional native entries beyond the previous 292-entry checkout:
  transferred extended Montgomery addition; conjugate and updated co-Z
  addition; Weierstrass and Montgomery signed recovery; three quotient relation
  residuals; Edwards w doubling and differential addition; a general
  Weierstrass chord formula; Hamburg's Figure 3 coupled ladder step.

The [paper import records](../research/paper-imports.json) give the exact source
and operation contract for each of those twelve entries. Derived/homogenized
formulas are labeled as such. These are not claims that the formulas are novel
or that every imported implementation is operation-count optimal.

## Equations and transformations

Every indexed model exposes its native general-Weierstrass forward/inverse maps
and a-invariants. On an ordinary affine chart those encode

```
a0*v^2 + a1*u*v + a3*v = u^3 + a2*u^2 + a4*u + a6
```

An a0-invariant other than one matters: Jacobi quartic's native map has a0=2.
The eight conditional inter-model maps, isogeny recipe and concrete target-curve
restrictions remain in the [atlas](../transformations/catalog.json). Model
availability is conditional on the field, roots, nonsingularity and torsion
structure; importing a family does not assert that every target curve admits it.

The generalized-Huff map uses `u=1/(x*y)`,
`v=(1+b*x*y)/(x*y^2)` and maps to
`v^2=u*(u+a)*(u+b)`. The inverse is
`x=v/(u*(u+b))`, `y=v/(u*(u+a))`. This is an independently derived open-chart
map, symbolically checked against the defining curve, with `a*b*(a-b) != 0`.
Excluded divisors require other charts; in particular the Huff neutral point
(0,0) is not in the forward affine chart. The model family is attributed to
[Wu–Feng](https://eprint.iacr.org/2010/390).

## Verification levels

The current index distinguishes:

| Status | Entries | Evidence |
| --- | ---: | --- |
| `finite_field_checked` | 351 | 200 recovered prime-field formulas, 135 binary formulas, twelve new entries and four earlier arithmetic additions |
| `inherited_not_rechecked` | 288 | Existing export formulas retained without a new independent arithmetic check |

The 200 recovered prime-field formulas passed 12,315 distinct checked cases over multiple
small prime fields. Their metadata and independently compiled three-operand
programs are both executed, and outputs are compared with independent
Weierstrass arithmetic. Inputs include equality constraints, parameter
specializations, projective rescalings and cached-coordinate invariants.

The 135 binary formulas passed 6,744 distinct checked cases over GF(16), GF(32) and GF(64),
using polynomial-basis extension fields and the general characteristic-2
Weierstrass law. Integer constants are reduced in characteristic 2, while
field-element integer encodings remain polynomial-basis values. Arithmetic
modulo 16 or 32 is **not** the field arithmetic used here.

The reports record every formula's case count, fields and chart/denominator
exclusions, with metadata and operation-code SHA-256 hashes:

- [Prime-field checks](../research/import-verification.json)
- [Binary checks](../research/binary-verification.json)
- [Paper-entry checks and contracts](../research/paper-imports.json)

The twelve paper/derived entries have independent group tests. The extended
Montgomery tests include neutral and 2-torsion points under complete-parameter
conditions. Hamburg's step is checked against affine P,Q,R=P+Q triples and
updated denominators. Co-Z tests check both outputs and their shared Z. Recovery
checks full points, and relation tests check nontrivial residuals rather than
only accepting zero. The extended linear coordinate maps are also checked
symbolically.

The native general-Weierstrass curve identity and inverse composition are
symbolically checked for all 13 prime-field model definitions, reducing
identities modulo the source curve ideal. Affine lifts and homogeneous weights
are checked across all 40 prime-field representations. These are open-chart
checks: they are not implementations of every exceptional projective chart.

Finite-field testing is evidence, not a symbolic proof for every field or a
proof of completeness. A skipped chart is counted, not silently treated as a
successful case. A wrong result on a tested valid chart fails the verifier.
`unified` is retained as a source tag and is never interpreted as `complete`.
No ECDLP, index-calculus or constant-time implementation speedup is established
by these tests.

## Joint-state and relation interfaces

The original operation types remain. New operation directories have explicit
contracts in the index:

| Operation | Input/output meaning |
| --- | --- |
| `conjugate` | Shared-Z P1,P2 produce P3=P1+P2 and P4=P1-P2 with Z3=Z4 |
| `updateadd` | Shared-Z P1,P2 produce P3=P1+P2 and rescaled P4=P1 with Z3=Z4 |
| `recovery` | Signed affine anchor P1, x(P2), x(P1+P2) produce full P4=P2 |
| `relation` | Quotient coordinates produce scalar residual R; not a point |
| `coupledstep` | A Hamburg joint ladder state produces its next state, not an arbitrary point addition |

Hamburg's entry is a conditional arithmetic step with explicit denominator
tracking. It is not a complete scalar-multiplication routine: scalar recoding,
swaps, setup/finish, exceptional-point avoidance and protocol behavior remain
separate. Likewise, invariant w loses additional torsion cosets and must not be
used as though it were a full-point coordinate.

Consumers that recognize only the original EFD operation types should opt in to
these directories using their contracts. Existing four-family consumers need
not load the new characteristic-2 namespace.

## Reproducibility

```
python -m pip install -r tests/requirements.txt
python -m unittest discover -s tests -v
python tools/verify_imports.py
python tools/verify_binary_imports.py
python tools/build_database.py
python tools/validate_database.py
python tools/search_formulas.py --field prime --source Hamburg --json
python tools/search_formulas.py --field binary --model shortw --verified
```

The import manifests retain exact upstream URLs and raw-source hashes:
[prime-field](../research/import-manifest.json) and
[binary](../research/binary-import-manifest.json). Products are made explicit,
non-identifier temporary names are renamed, equality notation is normalized,
and differential-addition outputs use index 4. Arithmetic programs are compiled
with an allowlisted AST; remote expressions are never passed to Python eval.
The index builder rejects stale arithmetic evidence if the file hashes change.

The live audit can be rerun without overwriting the historical gap snapshot:

```
python tools/audit_upstream_efd.py --output research/upstream-efd-current.json
```

## Remaining coverage

This establishes complete path coverage of both audited upstream genus-one EFD
sections, plus the stated additions. It does not establish a worldwide ranking
or exhaustive coverage of every research paper. The [17-paper backlog](../research/paper-gaps.json)
records which work is partially implemented and which remains a candidate.
Outstanding work includes the scaled optimal 2024 Montgomery laws, other Hamburg
states and complete ladders, LadderPrime's whole algorithm, further Edwards
invariants, higher summation-resultant recurrences, mu4/Lyness embeddings and
broader isogeny maps. Restricted-access papers remain identified as such.

Higher-genus Jacobians, supersingular characteristic-2 families outside the
ordinary-binary section and curves over other characteristics are not included
in the present coverage totals. Formula aliases, twists and parameter choices
must be deduplicated mathematically before making claims about the number of
independent schemes.
