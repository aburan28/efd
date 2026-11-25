# Additional arithmetic formulas and a rational 2-isogeny

These additions complement the [birational atlas](curve-transformations.md).
They are established identities or explicit algebraic derivations, not claims
of new cryptographic algorithms. Native arithmetic entries retain the existing
metadata/`.op3` pairing and output indexing.

## RFC 7748 Montgomery arithmetic

[RFC 7748 section 5](https://www.rfc-editor.org/rfc/rfc7748#section-5) specifies
`a24=(a-2)/4`, and a doubling step with

```
AA = (X+Z)^2
BB = (X-Z)^2
E = AA-BB
Xdouble = AA*BB
Zdouble = E*(AA+a24*E)
```

Existing entries use `a24=(a+2)/4` and `E*(BB+a24*E)`. These expressions are
equal with their respective constants, since `AA-BB=E`. Constants and expressions
must travel together. For Curve25519 the RFC constant is 121665; the older
convention gives 121666. No arithmetic speedup is claimed for switching forms.

Added native entries:

- [dbl-rfc7748](../montgom/xz/doubling/dbl-rfc7748): doubling.
- [ladd-rfc7748](../montgom/xz/ladder/ladd-rfc7748): homogeneous ladder with a
  projective difference `(X1:Z1)`; inputs P2, P3 give P4=2P2 and P5=P2+P3.
- [mladd-rfc7748](../montgom/xz/ladder/mladd-rfc7748): the RFC-style mixed case,
  with the explicit assumption `Z1=1`.

For the generic ladder, the addition output is

```
X5 = Z1*((X3-Z3)*(X2+Z2)+(X3+Z3)*(X2-Z2))^2
Z5 = X1*((X3-Z3)*(X2+Z2)-(X3+Z3)*(X2-Z2))^2
```

These are differential formulas, not arbitrary complete point addition.
Degenerate differences, including infinity or the rational 2-torsion point,
need separate handling; an output `(0:0)` is invalid. Protocol scalar processing,
conditional swaps, encoding rules, clamping and constant-time implementation
remain outside the formula files. Tests compare both metadata and operation code
against independent affine arithmetic on three prime fields, including nonunit
Montgomery `b` and projective rescalings.

## Short-Weierstrass XZ tripling

[tpl-division](../shortw/xz/tripling/tpl-division) adds the previously absent
`shortw/xz/tripling` operation using division polynomials. For
`y^2=x^3+a*x+b`, write

```
psi2 = 2*y
psi3 = 3*x^4+6*a*x^2+12*b*x-a^2
psi4 = 4*y*(x^6+5*a*x^4+20*b*x^3-5*a^2*x^2-4*a*b*x-8*b^2-a^3)
x(3P) = x - psi2*psi4/psi3^2
```

These are the standard multiplication-map division-polynomial identities (see
Silverman, *The Arithmetic of Elliptic Curves*, chapter III). Eliminating `y^2`
with the curve equation and homogenizing `x=X/Z` gives

```
H = 3*X^4+6*a*X^2*Z^2+12*b*X*Z^3-a^2*Z^4
G = X^3+a*X*Z^2+b*Z^3
F = X^6+5*a*X^4*Z^2+20*b*X^3*Z^3-5*a^2*X^2*Z^4
    -4*a*b*X*Z^5-(8*b^2+a^3)*Z^6
Xtriple = X*H^2-8*G*F
Ztriple = Z*H^2
```

Both output polynomials are homogeneous of degree 9. Use a nonsingular curve
in characteristic greater than 3. The formula also handles infinity and
nonzero rational 3-torsion as projective outputs; it does not recover the y
coordinate. Tests enumerate every quotient point on three small curves,
including infinity and torsion, with two projective scales, and compare against
threefold affine addition. It is a baseline with no claim of optimal operation
count; specialized tripling formulas may be faster.

## Rational 2-isogeny and its dual

The [separate isogeny record](../transformations/two-isogeny.json) gives the
classical degree-two map, a special case of Velu's formulas, on

```
E:  y^2 = x^3+A*x^2+B*x
E': y^2 = x^3-2*A*x^2+(A^2-4*B)*x
phi(x,y) = (x+A+B/x, y*(1-B/x^2))
```

Require `B*(A^2-4*B) != 0` and characteristic greater than 3. Direct substitution
proves the target identity. The kernel is `{O,(0,0)}`, both mapping to infinity.
For general short Weierstrass, translating a rational cubic root `r` to zero
produces `A=3*r`, `B=3*r^2+a`. Unlike Montgomery normalization, this translation
and isogeny do not require a square root of `B`.

Apply the same recipe to `E'`, then scale the resulting coordinates by `1/4`
and `1/8`, to obtain the dual

```
phihat(x,y) = ((x-2*A+(A^2-4*B)/x)/4,
              y*(1-(A^2-4*B)/x^2)/8)
phihat(phi(P)) = [2]P
```

The tests prove the target curve identity and dual composition symbolically
and check the homomorphism on finite-field points, including the kernel. On an
odd-prime-order subgroup, transport is injective and inverse(2) modulo the
subgroup order can undo the dual composition. It is not invertible on the full
even-order group. Curve25519 supports this recipe directly (`A=486662,B=1`);
P-256 and secp256k1 lack rational 2-torsion over their base fields.

The target has a general x coefficient `A^2-4*B` and is not automatically a
normalized Montgomery model. Reusing Montgomery formulas without checking a
normalizing square root would be incorrect. The record stays separate from
the repository's four native curve models for that reason.

## Verification

```
python -m pip install -r tests/requirements.txt
python -m unittest discover -s tests -v
```

`tests/test_additional_formulas.py` also checks that every `.op3` line contains
one binary operation and that metadata computations match the corresponding
operation code. These algebraic tests do not audit machine-code timing behavior
or establish an index-calculus speedup.
