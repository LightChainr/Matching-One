# P267 fixed-K density-clock orthogonalization

## Outcome

Exact density-clock and conditional-source removal eliminates about **68%** of
the raw external Euler response, but a highly resolved complex P4 response
remains.  The independent Target 1 and two-observer production blocks reproduce
that remainder at both sizes.

The retrospective protocol was committed as `fe26a8f` before this score, and
the scorer plus exact preflight was committed as `285229b` before numerical
reveal.  All roots, fixed-K conditional means, `beta_raw`, and `beta0` were
recomputed inside every delete-one replicate.

| independent block | N | raw P4 | fixed-K residual P4 | retained fraction |
|---|---:|---:|---:|---:|
| Target 1 | 325 | -22.0351 + 21.9598i | -7.07859 + 7.05511i | 0.321257 ± 0.000570 |
| Target 1 | 425 | -27.0494 - 18.1551i | -8.73512 - 5.84810i | 0.322680 ± 0.000575 |
| two-observer | 325 | -22.1014 + 21.8957i | -7.08584 + 7.04126i | 0.321090 ± 0.000536 |
| two-observer | 425 | -27.0261 - 18.1120i | -8.73217 - 5.85280i | 0.323115 ± 0.000543 |

Every residual complex vector rejects zero far beyond the frozen alpha=0.01;
the corresponding 2D chi-square statistics range from 184,417 to 223,841.
The underflow-safe log10 p-values are stored in `score.json`.

## Transfer and independent reproduction

| block | amplitude N425/N325 | phase (rad) |
|---|---:|---:|
| Target 1 | 1.05183 ± 0.00336 | 1.37369 ± 0.00513 |
| two-observer | 1.05233 ± 0.00316 | 1.37273 ± 0.00518 |

The two independent counter blocks are compatible for the fixed-K response:

- N325: chi-square=0.1169 on 2 coordinates, p=0.9432;
- N425: chi-square=0.00966 on 2 coordinates, p=0.9952.

Thus the density clock is a large component of the original bridge, but it is
not the whole bridge.  Approximately one third of the raw complex amplitude is
intrinsic to fluctuations within a fixed occupation layer and survives the
same-path `D0-beta0*S0` projection.

## Exact conditioning rule

For each occupation layer,

\[
\mu_N(k)=k-2N\frac{(k)_2}{(N)_2}+N\frac{(k)_4}{(N)_4}
             =E[O_{ext}\mid K=k].
\]

The scored source variables are `D0=JD-E[JD|K]` and
`S0=JS-E[JS|K]`, with `beta0=Re E[D0 conj(S0)]/E[|S0|^2]`.
The exact expansion and its simplified within-layer covariance agree to below
`9e-50` in every reported fit.  The complete delete-one covariance, including
raw/residual cross-covariance and retained fraction, is stored in the JSON.

## Boundary

This is a retrospective mechanism decomposition of two existing production
blocks.  It is not a new exponent fit or a continuum-field identification.
Only `O_ext` has the exact fixed-K mean used here.  No conditional mean or
density-clock score was invented for `O_far`.

