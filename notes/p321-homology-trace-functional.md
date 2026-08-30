# P321: the homology contrast is a two-closure trace, not one local central defect

## Result

At periodic connectivity widths 2, 3 and 4, there is a minimal exact
realization of the formal topological contrast

```text
Phi_Q(A) = Z_cross(A) - Q Z_trivial(A)
         = Tr(T A) - Q Tr(A)
         = Tr((T-QI) A).
```

Here `A` is a word in periodic FK connectivity joins, `T` is one-site
translation, and the two terms are respectively crossed-closure and ordinary
closure traces.  The internal matrix representative is therefore

```text
D_hom = T-QI.
```

The exact small-width certificate gives a sharp classification:

1. `D_hom` is not an ordinary central element from the first nondegenerate
   width 3;
2. `T` is an exact crossed-channel seam, but `T-QI` is not one crossed defect;
3. the contrast is naturally a difference of two closure traces, or a central
   sector label only after the two closure modules are adjoined as an external
   direct sum.

Thus the P321 missing object is genuinely module-resolved.  It cannot be
replaced by one local operator already living in the ordinary connectivity
module.

## Minimal periodic connectivity representation

Use the Catalan basis of circular noncrossing set partitions at width `w`.
The exact basis dimensions are

```text
w = 2,3,4
dim = 2,5,14.
```

Let `e_i` be the idempotent connectivity join of cyclic neighbours `i` and
`i+1`.  Translation acts by rotating the site labels.  Integer matrices give
the exact identities

```text
T^w = I,
T e_i T^-1 = e_(i+1),
T e_i - e_(i+1) T = 0.
```

The last line is the crossed-seam pull-through law.  It has zero matrix
residual at every site for all three widths.

For the homology contrast representative,

```text
[D_hom,T] = 0,
[D_hom,e_i] = [T,e_i],
D_hom e_i - e_(i+1)D_hom = Q(e_(i+1)-e_i).
```

The residuals are exact polynomials in `Q`, not floating-point norms:

| width | dimension | rank `[D_hom,e_i]` | squared norm | rank pull-through coefficient | squared norm coefficient |
|---:|---:|---:|---:|---:|---:|
| 2 | 2 | 0 | 0 | 0 | 0 |
| 3 | 5 | 2 | 6 | 2 | `6 Q^2` |
| 4 | 14 | 6 | 20 | 6 | `20 Q^2` |

Width 2 is not positive evidence for centrality: translation acts identically
on both connectivity states, so the seam is invisible.  Width 3 is the first
valid obstruction.

## Why the answer is a trace difference

The crossed component obeys the exact twisted trace law.  If

```text
sigma(A) = T A T^-1,
tr_T(A) = Tr(T A),
```

then

```text
tr_T(A B) = tr_T(sigma^-1(B) A).
```

The machine certificate checks this with zero integer residual on every pair
of generators.  Ordinary `Tr` is instead cyclic.  Their weighted difference
obeys neither single law:

- at width 3, `A=e0`, `B=T^2 e0` gives the ordinary cyclicity residual `-1`;
- at width 4, `A=e0`, `B=T e0` gives the same residual `-1`;
- with `A=B=e0`, the crossed-law residual of `Phi_Q` is exactly `-Q` at width
  3 and `-3Q` at width 4.

So `T-QI` is useful as the finite matrix which evaluates the difference, but
it must not be promoted to a transparent defect.  The two summands have
different closure laws.

There is a tautological central realization on

```text
M_cross direct_sum M_trivial
```

using the external label `diag(+I,-QI)`.  That construction does not recover a
hidden center of the ordinary module; it explicitly adjoins the additional
homology sector that P321 was trying to identify.

## Consequence for `F_t(tau)`

This certificate names the type of the missing object but does not determine
its modulus dependence.  A physical calculation still needs calibrated
generic-`Q` transfer weights and homology bookkeeping so that

```text
Z_cross = Z_2D,
Z_trivial = Z_0D
```

in the intended FK normalization.  The thermal one-point then comes from a
thermal insertion/derivative in each module before taking the trace
difference.  None of its spectral amplitudes, modular dependence or spinful
Ward data follows from the finite-width commutator identities.

Therefore the exact outcome strengthens, rather than removes, the obstruction
at `c2a5e2d`:

```text
F_t(tau)=<epsilon [1_(2D)-1_(0D)]>_tau
```

requires two sector-resolved traces.  The leading `E4` curve still cannot fix
the identity-dressed shape without this extra module data.

## Scientific card

- **Mechanism space changed:** excludes an ordinary central `D_hom` and a
  single transparent crossed defect in the minimal periodic connectivity
  module; retains a two-closure/direct-sum realization.
- **Not proved:** no identification of the formal crossed/ordinary closures
  with the fully weighted physical `Z_2D/Z_0D`, no value of `F_t(tau)`, and no
  tube-algebra classification.
- **Observer/sector/source/geometry:** exact homology trace functional;
  crossed versus trivial periodic closures; widths 2, 3 and 4; no stochastic
  source.
- **Dependency group:** exact algebraic continuation of P321 commit `c2a5e2d`;
  it reuses the noncrossing codec and no Monte Carlo archive.
- **Next lifting observation:** construct the smallest generic-`Q` affine-TL
  transfer word whose loop/homology bookkeeping independently calibrates both
  closures, then insert the thermal edge operator in each trace.  A full tube
  algebra is not required for that next gate.

## P370 proof-carrying interpretation

This is a small Level-0/Level-2 certificate in the sense of Issue 370:

- the identities and obstructions are exact integer or polynomial matrices;
- width 2 is explicitly tagged as a false-positive degeneracy;
- widths 3 and 4 provide portable rank and trace-law witnesses;
- the surviving mechanism is typed as an added module, not as a larger
  unconstrained latent fit.

The executable oracle is `scripts/p321_homology_trace_certificate.py`; the
frozen machine certificate is
`analysis/p321_homology_trace_certificate.json`.
