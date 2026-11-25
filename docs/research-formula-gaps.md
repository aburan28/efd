# Literature and upstream coverage gaps

**Status update, 2026-11-12:** all 200 upstream paths identified below are now
imported and independently finite-field checked. The complete ordinary-binary
EFD section and twelve further paper-derived/derived entries have also been
added. See [current coverage](coverage-and-verification.md). The dated inventory
below preserves the discovery rationale; paper records now identify partial
implementations and remaining work.

Checked 2026-11-09 against commit `6f1af53` of this export. This is a discovery
inventory, not a claim that every candidate is mathematically new or ready for
import. [Paper records](../research/paper-gaps.json) retain access status,
formula locations, overlap, field restrictions and target-curve applicability.
[The upstream snapshot](../research/upstream-efd-inventory.json) lists every
missing upstream formula path.

## The measurable gap

The [full large-characteristic EFD](https://www.hyperelliptic.org/EFD/g1p/index.html)
contains 471 formula links, 33 representations and ten model families. This
checkout contains 292 native formula metadata files, including the recent PR
additions. Comparing exact relative paths finds **200 absent upstream entries**:

| Omitted family | Equation / defining structure | Representations absent | Formula paths absent |
| --- | --- | ---: | ---: |
| Doubling-oriented Doche–Icart–Kohel | `y^2=x^3+a*x^2+16*a*x` | 1 | 13 |
| Tripling-oriented Doche–Icart–Kohel | `y^2=x^3+3*a*(x+1)^2` | 1 | 15 |
| Hessian | `x^3+y^3+1=3*d*x*y` | 2 | 25 |
| Jacobi intersection | `s^2+c^2=1`, `a*s^2+d^2=1` | 2 | 40 |
| Jacobi quartic | `y^2=x^4+2*a*x^2+1` | 6 | 90 |
| Twisted Hessian | `a*x^3+y^3+1=d*x*y` | 1 | 11 |
| Existing short Weierstrass family | Six old Chudnovsky addition paths | 0 | 6 |

Thus 194 absent entries are in omitted families; six are in an existing family.
The export README explicitly describes removing broken operation code, including
`shortw/projective-1/addition/add-1986-cc`. The six old Weierstrass entries require
individual review, rather than automatic restoration. Missing paths are not
necessarily missing mathematical identities. Conversely, identical paths do
not establish identical upstream and local content.

The snapshot has 271 paths in common and 21 local paths not linked by upstream,
which reconciles `471-271=200` and `292-271=21`. The 21 include local negations,
scaling and recent arithmetic additions. We did not download or execute the 200
missing formulas. The linked metadata and operation-code URL templates support
subsequent staged imports.

Reproduce the inventory (network access required):

```
python tools/audit_upstream_efd.py --output research/upstream-efd-inventory.json
```

The script follows the index's model and coordinate pages, deduplicates metadata
links, bounds concurrency at four requests and fails if a page cannot be read.
Counts come from discovered links, not from subtracting potentially stale
homepage totals. The frozen JSON records this particular audit date.

## Missing candidates from research papers

The strongest import priorities are extended Montgomery embeddings, coupled
short-Weierstrass ladders, shared-denominator operations and quotient/recovery
relations. They add capabilities absent from the current formula interfaces.
Paper-reported costs below are not measurements from this repository.

| Paper | Concrete candidate to extract | Overlap and applicability |
| --- | --- | --- |
| [Farashahi–Fadavi–Sabbaghian, 2024](https://ches.iacr.org/2024/papers-issue-4/4_111.pdf) | Extended Montgomery coordinates; linear Edwards maps, Eqs. (33)–(34); scaled complete laws (44)–(45) | `montgom` currently has only XZ. Table 3 and the Appendix identify the conditional `8M+1D` law for Curve25519. This does not apply to P-256 as a Montgomery model. |
| [Hamburg, 2020](https://www.shiftleft.org/papers/ladder/ladder-tches.pdf) | Figures 3–5: two coupled Montgomery-ladder states and a Joye-ladder state on short Weierstrass curves | Direct P-256 candidate. The first ladder step is reported as `8M+3S+7A`; setup, recovery and completeness have separate costs. [Author code](https://github.com/bitwiseshiftleft/ladder_formulas) is available. |
| [Hışıl–Pala, 2026, LadderPrime](https://dergipark.org.tr/en/download/article-file/5622508) | Theorems 1–4 and Algorithm 1: general XZ doubling and additive/multiplicative differential formulas | Direct generic-Weierstrass/Kummer candidate. Audit the complete algorithm's exception handling; individual differential formulas have restricted inputs. |
| [Goundar–Joye–Miyaji–Rivain–Venelli, 2011](https://www.matthieurivain.com/files/jcen11b.pdf) | ZADDU, ZADDC, DBLU/TPLU, fused steps and XY-only states | Meloni 2007 `zadd` already exists. The missing capabilities are updated shared-Z outputs, simultaneous sum/difference and fused state transitions. Relevant to P-256. |
| [Chung–Costello–Smith, 2015](https://arxiv.org/abs/1510.03174) | Sections 5.1–5.3: signed recovery for Montgomery, Edwards and Weierstrass; two-dimensional pseudomultiplication | Many basic operations already exist; recovery is the interface gap. Covers both target curve families with the appropriate model and signed anchor. |
| [Renes–Smith, 2017, qDSA](https://joostrenes.nl/publications/qdsa-eprint.pdf) | Section 3.2, Proposition 3: biquadratic sum/difference relation; Algorithm 2: quotient relation checking | Direct Curve25519 candidate. An unordered quotient relation, not an ordinary addition formula. The relation-profile tool covers related algebra but does not implement this operation. |
| [Hosseini–Farashahi, differential Edwards arithmetic](https://arxiv.org/pdf/2606.20831) | Section 4: alternative invariant w-coordinates, differential arithmetic and recovery | The 2017 conference paper has a 2026 arXiv deposit; deduplicate versions. Includes `w=d*x^2*y^2`, with extra torsion-coset ambiguity. Existing Edwards square-d YZ formulas do not cover all these variants. |
| [Cheng–Kodera–Miyaji, 2019](https://www.jstage.jst.go.jp/article/transfun/E102.A/9/E102.A_1061/_article) | Model-specific summation-polynomial comparisons | Closely matches the coordinate/index-calculus question. Publisher abstract checked; full text requires access. No exact formulas, field-regime conclusion or timing claims imported from the abstract. |
| [Bianco–Gorla, 2015 preprint](https://arxiv.org/pdf/1507.08555) | Definition 9: Edwards summation polynomial and resultant recursion; Theorem 10; Remark 14 | Supplies explicit provenance and higher-order relations beyond a transported f3. Trace-zero compression experiments concern extension fields, not a demonstrated prime-field ECDLP improvement. |
| [Hışıl–Renes, 2019](https://joostrenes.nl/publications/kummerlines.pdf) | Propositions 3.1–3.2, Theorem 3.4: canonical/squared Kummer, Montgomery and Edwards map network | Some older theta/YZ arithmetic overlaps EFD. Full rational 2-torsion is required; Curve25519 and P-256 themselves do not meet it. Sharing a base prime is insufficient. |
| [Kohel, geometry of efficient arithmetic](https://arxiv.org/pdf/1601.03665) | Projective-embedding/addition-law search; Section 6, Theorem 6.11: special doubling families | Broader search framework. General embeddings can apply to P-256; the specialized 2-/3-torsion optimizations require extra structure. |
| [Kohel, twisted mu4-normal form](https://arxiv.org/abs/2012.10799) | Odd-characteristic linear Edwards isomorphism; alternative intersections of quadrics | A missing representation family. The headline operation costs in the abstract concern binary fields and must not be applied to P-256. |
| [Joye–Tibouchi–Vergnaud, 2010](https://eprint.iacr.org/2010/383) | Huff models, maps, addition/doubling | Primary abstract checked. Standard Huff needs a `Z/4 x Z/2` subgroup; neither target curve has the required full rational 2-torsion. |
| [Wu–Feng, 2010](https://eprint.iacr.org/2010/390) | General Huff `x*(a*y^2-1)=y*(b*x^2-1)` and projective operations | Primary abstract checked. Full rational 2-torsion is again required. Broader atlas entry, not a base-field rewrite of either target. |
| [Hone, 2020](https://arxiv.org/pdf/2002.03811) | Lyness/QRT biquadratic curves in `P1 x P1`; Sections 2–4: maps, translation and doubling | A genuinely different embedding/search direction. No explicit finite-field conversion for the target curves was established here. Fixed-point translation is not arbitrary addition. |
| [Moody, 2011](https://www.nist.gov/publications/division-polynomials-jacobi-quartic-curves) | Jacobi-quartic division-polynomial recurrences and general multiples | NIST record checked. The family is omitted from this export; exact recurrence import still requires full-text review. |
| [Amadori–Pintore–Sala, 2018](https://pure.tue.nl/ws/files/91158650/609.pdf) | Prime-field point decomposition and modified index-calculus construction | A relation-generation baseline, rather than a new coordinate model. Useful alongside representation experiments. |

## An example of a missing representation

The 2024 extended-Montgomery paper gives a concrete addition to the map atlas.
Using its coordinate order, extended Edwards is

```
a*X^2 + Y^2 = Z^2 + d*T^2
X*Y = Z*T
```

and extended Montgomery is

```
B*S*V = U^2 + A*U*W + W^2
S*U = V*W
A = 2*(a+d)/(a-d), B = 4/(a-d)
```

Its linear maps are

```
(U,V,S,W) = (X+T, Z+Y, Z-Y, X-T)
(X,Y,T,Z) = (U+W, V-S, U-W, V+S)
```

Composing the maps multiplies all coordinates by 2, which is the same
projective point in odd characteristic. This is a different projective
embedding from Montgomery XZ: it retains the full point. The affine chart
has `u=U/W`, `v=V/W`. The paper also derives scaled variants for cheaper laws;
those variants and their character conditions must be imported together.
These transcribed identities are research candidates, not native validated
formula entries in this change.

Another important candidate is `w=d*x^2*y^2` from the differential-Edwards
paper. It removes more information than ordinary y: it identifies additional
torsion cosets. Any relation derived in w needs lift/subgroup validation and
cofactor accounting. Comparing only the number of symbols or the apparent
degree would miss that cost.

## How to organize a larger database

Count model families, coordinate representations, arithmetic formulas and
parameter specializations separately. Repeated transformations and rescalings
create unlimited presentations; they do not automatically yield independent
algorithms. A useful taxonomy is:

- **Equations:** long/short Weierstrass, Montgomery, Edwards/twisted Edwards,
  Hessian/twisted Hessian, Jacobi quartic/intersection, both DIK families,
  Legendre/full-2-torsion forms, Huff/general Huff, mu4/twisted mu4, and
  Lyness/QRT biquadratic models.
- **Embeddings:** plane cubics, weighted projective presentations, intersections
  of two quadrics in P3, product-projective `(2,2)` curves, and higher-degree
  elliptic normal embeddings. Kohel's framework supplies a systematic search
  through equations modulo the embedding's ideal.
- **Full-point representations:** affine, homogeneous, Jacobian, modified and
  Chudnovsky-cached Jacobian, XYZZ, inverted, extended Edwards, extended
  Montgomery and family-specific cached quartic/intersection coordinates.
- **Joint states:** shared-Z/co-Z, shared powers, mixed affine/projective,
  updated-addition, conjugate addition and coupled ladder states. These require
  multi-output interfaces and state invariants; merely adding another
  `variables` file does not capture them.
- **Quotients/relations:** XZ, YZ, squared YZ, theta/Kummer lines, invariant
  w-functions, biquadratic relation checks and summation-polynomial variables.
  Record the map degree, fibers, lost signs/torsion and recovery data.
- **Different field regimes:** binary Edwards, Lopez–Dahab, lambda and related
  binary-field schemes belong in a separate characteristic-2 catalog. They are
  not additional prime-field choices just because their coordinates look similar.

For P-256, prioritize arbitrary projective embeddings, co-Z/coupled ladders and
Weierstrass quotient relations. Rational-2-torsion forms are excluded over the
base field; Hessian models need rational 3-torsion. For Curve25519, prioritize
extended Montgomery and the alternative Edwards invariants. Its rational
2-torsion does not imply *full* rational 2-torsion. The tripling-oriented DIK
family involves rational 3-isogeny structure; do not equate a rational isogeny
kernel with a rational point of order 3 without checking Frobenius action.

For index calculus, use the summation-polynomial papers to compare the same
factor base and subgroup across presentations, with denominator saturation,
quotient-fiber checks and total decomposition cost. Fast addition alone does
not establish a better relation solver.

## Import order

First recover the 194 omitted-family entries with source provenance and the
existing export's normalization/verification rules. Separately review the six
old Weierstrass paths. Then add the 2024 extended model and joint-state/recovery
interfaces before importing their arithmetic. Keep every candidate tagged with
its source location, characteristic, parameter conditions, domain, completeness
meaning, operation-count convention and independent verification status.

The 17-paper catalog is a starting bibliography. It is not a claim of exhaustive
coverage of all published coordinate systems, nor a claim of 200 new schemes.
