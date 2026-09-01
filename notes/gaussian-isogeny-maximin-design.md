# Gaussian isogeny design: the four-edge maximin campaign

This note freezes a concrete experiment, not a general research agenda.  The
exact generator is `scripts/optimize_gaussian_isogeny_design.py`; the complete
targets and budgets are in
`predictions/gaussian_isogeny_maximin_design_20260828.yaml`.

## Result

Enumerating every primitive Gaussian multiplier of norm at most 65 from the
high-signal `N=65,85` parents leaves the existing norm-5 children at the top of
the balanced angular/radial information-per-cost ranking:

| rank | edge | child | main fingerprint |
|---:|:---|---:|:---|
| 1 | `65 --(1+2i)--> 325` | `(17,6)/(18,1)` | H4 negative, H12 positive |
| 2 | `85 --(1-2i)--> 425` | `(16,13)/(19,8)` | H4 negative, H12 positive |
| 3 | `65 --(1-3i)--> 650` | `(23,11)/(19,17)` | weaker norm-10 backup |
| 4 | `85 --(1+3i)--> 850` | `(29,3)/(27,11)` | weaker norm-10 backup |

The norm-5 exact angular factors are

\[
r_4=-\frac{14}{25},\qquad
r_8=-\frac{1054}{625},\qquad
r_{12}=+\frac{23506}{15625}.
\]

Thus H4 and H12 predict opposite child signs.  Larger primitive multipliers do
create attractive near-null harmonic edges—for example one norm-13 edge at
`N=1105` nearly nulls H12—but the child signal and site cost make them inferior
as the first experiment.

Norm 5 alone is not the best radial-exponent experiment.  The strongest small
campaign combines it with cheap norm-2 transfers.  On a 100M-replica grid and
a 737-billion-site-update budget, the frozen maximin allocation is:

| child | replicas | job |
|---:|---:|:---|
| 130 | 2.0B | radial leverage from N=65 |
| 170 | 0.6B | second-parent radial leverage from N=85 |
| 325 | 0.5B | H4/H12 sign split and S-prime transfer |
| 425 | 0.5B | independent H4/H12 sign split and S-prime transfer |

With the declared planning SE, the weakest pure-model separation is
`chi2=16.53` (about `4.07 sigma`) between radial exponents `4/3` and `9/8`.
The H4-versus-H12 separation is `chi2=31.38` (about `5.60 sigma`).  The same
threshold-rank histograms also score the frozen normalized `P4[S']` ratios:

| edge | pure | q=2 correction | Jordan log |
|:---|---:|---:|---:|
| 65 -> 325 | 0.1337481 | 0.2161628 | 0.2539065 |
| 85 -> 425 | 0.1337481 | 0.1871016 | 0.2382568 |

This is why the recommended compute is asymmetric: the small `N=130` child
buys radial information cheaply, while the two norm-5 children buy angular and
Jordan/correction information that norm 2 cannot supply.

## The norm-4 caveat is an engine boundary, not an optimality theorem

Multiplication by `2` has Gaussian norm four and gives the exact semigroup edge

\[
T_4=T_2^2
\]

up to the square-lattice-trivial unit `i`.  It was absent from the primitive
enumeration because multiplying a primitive parent by `2` produces content-2
periods:

- `N=65 -> 260`: `(16,2)/(14,8)`, Smith invariants `(2,130)`;
- `N=85 -> 340`: `(18,4)/(14,12)`, Smith invariants `(2,170)`.

The quotient groups are `Z/2 x Z/130` and `Z/2 x Z/170`, not cyclic `Z/N`.
The current C++ engines use a single cyclic label and therefore cannot run
these children.  Exact angular ratios are all `+1`, while the three radial
ratios are

\[
4^{-13/8}=0.1051121,\quad
4^{-4/3}=0.1574901,\quad
4^{-9/8}=0.2102241.
\]

Norm 4 is scientifically useful.  It creates the cheapest direct scale-
curvature triple `N -> 2N -> 4N`, checks the semigroup sign recovery, probes
whether cyclic and non-cyclic quotient groups share the same continuum
amplitude, and separates the frozen S-prime mechanisms.  For the N=65 lineage,
the pure/q2/Jordan normalized S-prime ratios are
`0.1767767/0.2788974/0.3135726`.

It is nevertheless not the best replacement for the primitive four-edge
campaign under the current fixed-p variance model.  Reoptimizing the same
budget with norm-4 plus norm-5 rows gives maximin `chi2=7.55`, versus `16.53`
for norm-2 plus norm-5.  The recommendation is therefore precise:

1. run/freeze the primitive norm-5 campaign without waiting;
2. implement a minimal Smith-coordinate general-period production backend in
   parallel;
3. use `N=260,340` as its first production target for radial curvature and
   backend universality, not as a substitute for H4/H12 discrimination.

The backend is worth implementing because it opens a genuinely new quotient-
topology axis.  It is not worth putting on the critical path of `N=325,425`.
