# Post-H4 residual: odd thermal spin4 x even identity spin4 as the first nonlinear mechanism

Date: 2026-09-14

Status: correction/addendum to `post-h4-residual-hierarchy-20260914.md`.  The earlier note placed a linear identity-family spin-eight candidate too high.  Once matching parity is enforced, the more natural first mechanism is the **second-order product of the already observed dual-odd thermal spin-four field with the already observed dual-even/common spin-four lattice anisotropy**.

This note changes candidate ordering; it does not erase the earlier exploratory calculation.

## 1. Two spin-four directions are already present

The fixed-width safe-transfer hierarchy contains two distinct finite-size structures.

### Common / dual-even spin four

The average magnetic excitation has a large correction consistent with

```text
x_even = 4,
y_even = 2-x_even = -2.
```

On a square lattice the natural lattice-anisotropy realization is a spin-four identity-family correction.  Its important feature for the matching root is not the exact field name but that it is large in the common/average channel and strongly cancels from the primal/matching difference.

### Difference / dual-odd spin four

The charge-sector difference has

```text
x_odd = 21/4,
y_odd = 2-x_odd = -13/4,
```

with strong `H4(theta)=cos(4theta)` evidence.  The ordinary thermal-Q4 Ward identity explains why its leading critical projection is thermal-tangent.

Thus the effective finite-size theory already contains

```text
v4^+ H4 ell^-2       (matching-even/common),
u4^- H4 ell^-13/4    (matching-odd/difference).
```

## 2. Matching parity changes the post-H4 candidate ordering

A linear identity-family spin-eight field is naturally matching-even unless a separate microscopic odd coupling is demonstrated.  Therefore it should not be the default first source of an odd charge-root residual.

By contrast, the product

```text
u4^- * v4^+
```

is automatically matching-odd:

```text
odd x even = odd.
```

It is therefore allowed in the charge difference without introducing any new primary/module.

## 3. Angular structure is locked: H4 squared

The mixed second-order term carries

```text
H4(theta)^2
 = [1+H8(theta)]/2.
```

So it generates a scalar and an H8 component with the **same field-level coefficient** in the standard `H0=1`, `H8=cos(8theta)` basis.

This is more predictive than a generic `S0+C8 S8` residual: the simplest mixed mechanism does not have two free angular amplitudes.

For a two-angle H4 projector with `h_i=H4(theta_i)`, the residual of `H4^2` is exactly

```text
P_perp4[H4^2]
 = (h1 h2^2-h2 h1^2)/(h1-h2)
 = -h1 h2.
```

Using

```text
C8=-(1+2h1h2),
```

this is equivalently

```text
P_perp4[H4^2] = (1+C8)/2.
```

Thus every existing same-circle pair has a parameter-free geometry factor for this mechanism.

## 4. Root exponent is exactly six at the scaling level

A second-order finite-size term built from the two couplings has RG exponent

```text
y_mix = y_odd+y_even
      = -13/4-2
      = -21/4.
```

An excitation-energy correction carries one extra inverse length:

```text
Delta I_mix
 ~ u4^- v4^+ H4^2 ell^(y_mix-1)
 = H4^2 ell^-25/4.
```

The thermal derivative which moves the charge root scales as

```text
partial_p Delta I ~ ell^-1/4.
```

Therefore the pseudo-critical root receives

```text
boxed:
delta p_mix
 ~ H4(theta)^2 ell^-6.
```

This provides a structural origin for the familiar next even power in the raw axial root sequence without introducing a new single field of root exponent six.

## 5. Why the H4 projector exposes this term

The same-circle H4 projector removes every correction proportional to `H4(theta)`, including a possible linear spin-four tower

```text
H4[a4 ell^-4+a6 ell^-6+...].
```

But it does **not** remove `H4^2`.

Therefore after leading angular improvement, a natural first surviving nonlinear contribution is precisely

```text
(-h1 h2) C_mix ell^-6.
```

This makes the projected residual scientifically useful rather than merely an estimator error: it directly probes nonlinear mixing between the leading odd and common even lattice directions.

## 6. Existing residuals do not yet reduce to one pure H4-squared term

For the current two-angle projectors, `h1` and `h2` generally have opposite signs, so

```text
-h1 h2 >0.
```

A **single** mixed coefficient with negligible competitors would therefore give the same residual sign across these pairs at comparable large `ell`.

The observed projected residuals include both signs:

```text
N65:   positive,
N85:   near zero positive,
ell=10: negative.
```

Hence the pure mixed term is not by itself a complete fit at current sizes.  This is valuable: it means the N1105 same-circle decomposition should not be reduced to one locked coefficient in advance.

Possible reasons include

```text
additional scalar odd channel,
independent H8 odd channel,
H12/higher harmonic,
higher-order finite-size contamination,
logarithmic/normal spin4 remainder.
```

The mixed mechanism remains structurally preferred as the **first nonlinear term**, not as an already sufficient numerical model.

## 7. N1105 gains a sharp mechanism test

The four same-circle N1105 roots can be decomposed as

```text
P0 + P4 H4 + P8 H8 (+ P12 H12).
```

The simplest odd-Q4 x even-spin4 mechanism predicts, for the post-H4 second-order piece,

```text
P0_mix = P8_mix
```

because both come from `(1+H8)/2`.

Thus define the same-ell mechanism ratio

```text
R_mix = P8/P0
```

only after the dominant P4 piece is separated.

Interpretation:

```text
R_mix ~ 1, P12 small:
    strong support for nonlinear odd/even spin4 mixing;

P8 nonzero, P0 small or R_mix far from 1:
    independent H8 channel required;

P0 dominates:
    separate angle-independent odd scalar/log channel;

P12 substantial:
    H0/H4/H8 truncation insufficient.
```

A single N cannot determine the root exponent six.  N1105 is an angular/parity mechanism test only.

## 8. Relation to the discarded simple coordinate-squared explanation

An earlier exploratory estimate considered the ordinary analytic nonlinearity produced by squaring the leading odd thermal-coordinate shift itself.

Symmetry makes the conceptual distinction clearer:

```text
(u4^-)^2: even under matching,
u4^- v4^+: odd under matching.
```

The charge difference is odd, so the first expression should not be the default source of the residual even before noticing that its calibrated numerical size is tiny.

The physically relevant second-order candidate is the **odd x even** mixed term.

## 9. Revised candidate ordering after H4 projection

Current ordering:

```text
1. odd thermal-Q4 x common even spin4 nonlinear mixing
   -> H0+H8 locked structure, root ell^-6;

2. additional independent dual-odd scalar/H8/log channels;

3. H12/higher D4 harmonics;

4. linear identity-family spin8 x=8
   only after a nonzero matching-odd coupling is independently justified.
```

A thermal-family spin-eight field remains a possible independent H8 channel, but its root exponent would be eight and its amplitude need not satisfy the locked H0=H8 relation.

## 10. Research consequence

The leading and next correction mechanisms may now be organized without inventing a new field for every power:

```text
linear odd thermal Q4 spin4
    -> root ell^-4 H4,

nonlinear odd-Q4 x even-identity-spin4
    -> root ell^-6 (H0+H8),

then genuinely new irreps/modules only if the locked second-order prediction fails.
```

This is a much more economical hypothesis than assigning the raw `4,6,...` ladder to unrelated operators.  It is also directly falsifiable by one same-circle multi-angle calculation.

## 11. Claim boundary

- The RG exponent arithmetic, matching-parity product and identity `H4^2=(1+H8)/2` are exact once the two input scaling fields are granted.
- Identifying the common `x≈4` correction specifically with a dual-even spin-four identity-family lattice field remains a scaling/CFT interpretation, albeit a standard one.
- The equality `P0_mix=P8_mix` is the leading second-order prediction; current two-angle residuals do not yet satisfy a one-term model cleanly.
- This note supersedes the earlier ordering that placed linear identity-spin8 first; it does not delete that field as a possible subleading contribution.
