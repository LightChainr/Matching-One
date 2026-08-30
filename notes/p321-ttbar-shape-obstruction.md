# P321: the leading E4 curve does not fix the identity-dressed correction

## Result

The P321 leading aspect-ratio result may be written conditionally as

```text
F_u(tau) = kappa E4(tau) F_t(tau).
```

This relation does **not** by itself produce a parameter-free prediction for
the next root coefficient `D_width/C_width` under an ordinary `T Tbar`
identity-family perturbation.  The missing object is not another lattice
amplitude.  It is the modulus dependence of the homology-resolved thermal
one-point `F_t(tau)`.

## Exact first-order operator

Expand the real-analytic torus one-point deformation in Cardy's Theorem 2
(arXiv:2201.00478).  For a Fourier mode

```text
exp(-2 pi x y + 2 pi i p theta)
```

of a scalar one-point function with scaling dimension `k`, the normalized
first-order multiplier is

```text
4 pi^2 (x^2-p^2) y^2 - 2 pi (1+k) x y.
```

Equivalently,

```text
L_k = y^2 (partial_y^2 + partial_theta^2) + (1+k)y partial_y.
```

This is the one-point analogue of the torus differential structure also seen
in the perturbative Ward calculation of He--Sun (arXiv:2004.07486).  It is
already enough to decide whether a cancellation follows from a product ratio.

For scalarized weights `kf,kg`, direct expansion gives

```text
L_(kf+kg)(fg) - g L_kf(f) - f L_kg(g)
 = 2 y^2 grad(f).grad(g)
   + kg y g partial_y(f) + kf y f partial_y(g).
```

Therefore

```text
L_(kf+kg)(fg)/(fg) - L_kf(f)/f
```

contains `grad log(f)`.  Substituting `f=F_t` and `g=E4` shows that the
subleading numerator/denominator ratio retains derivatives of `F_t`.  The
unknown leading lattice scalar coupling may be common, but the modulus shape
does not reduce to an E4 derivative alone.

## Stronger boundary from spin

`Q4 epsilon` is spin four.  Cardy's displayed real-analytic theorem applies to
a scalar one-point function.  A correct spinful formula requires modular
biweights/covariant derivatives.  The calculation above deliberately gives
the proposed cancellation its most favorable scalarized reading and it still
does not cancel.  A spinful treatment may specify the missing terms, but it
cannot be inferred from the leading ratio alone.

## Consequence for the research program

The next exact continuum object is now sharply named:

```text
F_t(tau) = <epsilon [1_(2D)-1_(0D)]>_tau
```

at `Q=1`, or its generic-Q defect/tube-algebra precursor.  The sphere--torus
connectivity-map or affine-TL defect program proposed in Issue 321 is not an
optional elegance: it supplies exactly the logarithmic derivative missing
from the identity-dressing prediction.

The deeper equal-area rectangle stream remains useful.  It can estimate a
stable empirical `D/C` curve and score any independently supplied `F_t`
model.  It should not be used to invent an arbitrary post-reveal E4-only
correction function.

## Claim boundary

This closes only the claimed automatic cancellation.  It does not reject
ordinary identity dressing, derive the spinful deformation, prove defect
transparency, or identify a logarithmic partner.  It converts those broad
alternatives into one concrete missing function.
