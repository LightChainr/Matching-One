# N100: a parameter-free three-modulus falsifier of affine-E4 coupling

Geometry input: `b9e4ea19bc585cbed18ec6ba1d13e85f2b5accc7`.
No Monte Carlo, archived-data refit, PR or comment is part of this result.

## The new observable

Let `Y2, Y4, Ys` be the four-vectors `P4(A_top,E_top,C,W)` on the N100
moduli `2i, 4i, 1/2+i`, respectively. Use the same fixed probability
`p_ref=.59274605079`, the exact first-minus-second cos4 normalization, and
the Smith pair `(1,100)->(5,20)` at every modulus. The six explicit matrices
and nonzero normalizations are fixed by `b9e4ea1`; this is not an arbitrary
continuous tau scan.

For each coordinate j allow completely unknown constants `a_j,b_j`, independent
of those for the other coordinates, and declare only

`Y_j(tau)=a_j+b_j*g(tau)`

at these three finite geometries. With `g(tau)=Im(tau)^2 E4(tau)/E4(i)`, the
following **four-vector zero is exact under this coupling hypothesis**:

`120 sqrt(2) Y2 + (47-60 sqrt(2)) Y4 - (47+60 sqrt(2)) Ys = 0`.     (E)

It removes every coordinate's additive offset and multiplicative amplitude.
It uses no N50 observed value, no inter-area exponent and no transfer scale.
It does not require the A/E/C/W vectors to share a ray.

## Two concrete adversaries give different zeros

For affine **height-only E4**, replace E4(tau) by E4(i Im(tau)). The shear is
then treated as an ordinary square of height one. Its exact null is

`(75+60 sqrt(2)) Y2 - 28 Y4 - (47+60 sqrt(2)) Ys = 0`.             (H)

For affine **height squared**, `Y_j=a_j+b_j Im(tau)^2`, the null is

`5 Y2 - Y4 - 4 Ys = 0`.                                         (Q)

These are fixed hypotheses, not three names for one fit. The equivalent
signed secant predictions are:

| hypothesis | r=(Ys-Y2)/(Y4-Y2) | numerical value |
|---|---|---:|
| affine-E4 | `(47-60 sqrt(2))/(47+60 sqrt(2))` | -0.287083852577789 |
| affine height-only E4 | `-28/(47+60 sqrt(2))` | -0.212358001359808 |
| affine y squared | `-1/4` | -0.250000000000000 |

**Score the linear zero, not an empirical noisy ratio.** With shear coefficient
one this is `R_f=Ys-Y2-r_f*(Y4-Y2)`. If model f is true, the expectation of
the null for another model h is exactly

`E[R_h]=(r_f-r_h) E[Y4-Y2]`.

All three slopes differ exactly. Their affine shape planes intersect only
in the constant three-shape response. Thus, if even one coordinate has a
nonzero 4i-minus-2i span, exact population data distinguish all three models.
An unresolved span instead gives an underpowered experiment, not evidence
for any mechanism.

The normalized separation is .074726 of that span between E and H, and
.037084 between E and Q. As a rough single-coordinate 3-SE planning guide,
the latter needs residual SE below .012361 times the span. This is not a
measured variance or an acquisition authorization. A joint four-coordinate
score can use the full covariance; it cannot multiply the four coordinates
as independent evidence.

## Why the coefficients are exact

Normalize `T=theta3(i)^4`; by the square elliptic symmetry
`theta4(i)^4=theta2(i)^4=T/2`, hence `E4(i)=3T^2/4`.
The elementary theta duplication identities imply

`theta3(2i)^4/T=(3+2 sqrt(2))/8 = u`,
`theta4(2i)^4/T=1/sqrt(2) = v`,

and

`E4(2 tau)=[theta3(tau)^8+14 theta3(tau)^4 theta4(tau)^4+theta4(tau)^8]/16`.

Applying this first at i and then at 2i yields

`g(2i)=11/4`,
`g(4i)=(4/3)(u^2+14uv+v^2)=(91+60 sqrt(2))/16`.

The index-two coset identity already derived in `b9e4ea1` is
`g(4i)+g(1/2+i)=91/8`, so

`g(1/2+i)=(91-60 sqrt(2))/16`.

Taking the cross product of `(1,1,1)` and each exact shape vector gives
(E), (H), and (Q). The new script performs these operations in Q(sqrt(2));
all offset and shape residuals are exactly zero. An independent q-series
evaluation is only a numerical cross-check, not the proof.

## Coupling boundary and executable covariance readout

The modular-form identities are unconditional mathematical statements.
Applying their affine span to finite `P4(A,E,C,W)` is the explicit physical
hypothesis. It requires each coordinate's offset and coupling to remain
constant across the three declared N100 geometries after the existing spin4
normalization. Geometry-specific counterterms, additional modulus functions,
or nonseparable finite-width effects can violate it. Such a violation rejects
this coupling law, not E4 as a function or a uniquely named continuum field.

Stack the three four-vectors shape-major into a 12-vector. The machine file
`predictions/p267-n100-three-modulus-null.json` supplies each exact weight,
decimal projection `L`, normalized slope and pairwise separation. The
companion `project_joint(mean,covariance,model)` returns

`R=L mean`, `Cov(R)=L Sigma L^T`.

Full same-stream cross-shape and cross-field covariance is required whenever
shared random input is used. A future joint Wald score has four residual
coordinates when the projected covariance has full rank; singularity is not
fixed by discarding correlations or pretending to gain independent samples.

Compared with the source-frozen two-cell shear/2i proposal, this three-cell
test costs one additional shape pair but removes the old-source uncertainty,
unknown new-area scale and cross-area shape-transfer assumption. This is a
scientific tradeoff, not a default request for extra production.

## Scientific card

- Changed mechanism space: three feasible same-area moduli make affine-E4,
  height-only E4 and affine y squared mutually identifiable except on the
  exactly flat response; a third shape can falsify the entire affine law.
- Not proved: that any law is physically selected, any target variance is
  affordable, or the observed global H4 is an E4/thermal/Jordan field.
- Observer/sector/source/geometry: four normalized rank/clock contrasts,
  Bernoulli square-site at fixed p, N100 cyclic/noncyclic pairs, same O.
- Dependency: exact theory only; future three-shape covariance is one block.
- Next discriminator: the frozen linear residual vector and its full
  covariance on the six declared period matrices, not a free shape refit.

Reproduction:

```
python3 scripts/p267_three_modulus_null.py --json predictions/p267-n100-three-modulus-null.json
python3 -m unittest discover -s tests -p test_p267_three_modulus_null.py -v
```

Four focused tests pass. No simulated samples were produced.
