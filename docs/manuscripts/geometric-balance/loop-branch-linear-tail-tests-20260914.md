# Linear-tail tests implied by the loop--branch variational candidate

2026-09-14.  Consequences of the **candidate** #758 rate function, not a proof that the actual SITE Palm span obeys that LDP.  The purpose is to turn the variational proposal into sharp falsification tests that can reuse existing rare-component machinery.

## 1. Candidate rate and its saturation point

The proposed fixed-subcritical NN component-Palm rate is

\[
I_p(A)=\min_{0\le r\le A}
\{\tau_p(1,2r)-\kappa+\kappa(A-r)\},                           \tag{1.1}
\]

with `kappa=tau_p(1,0)=tau_p(0,1)`.

As shown in `structural-consequences-20260914.md`, strict convexity of the directional norm implies a unique `r_*(p)>0`, independent of `A`, and

\[
I_p(A)=
\begin{cases}
\tau_p(1,2A)-\kappa,&A\le r_*,\\
\kappa A+c_*,&A\ge r_*.
\end{cases}                                                     \tag{1.2}
\]

Square symmetry gives the deterministic bound

\[
\boxed{0<r_*(p)<1/2.}                                         \tag{1.3}
\]

Therefore **every threshold `A>=1/2` is already in the saturated linear-branch regime of the candidate.**

## 2. Tail ratios for any two macroscopic thresholds beyond 1/2

If the actual Palm span obeys the candidate LDP

\[
-\frac1w\log P_{Palm}(L\ge\lceil Aw\rceil)\to I_p(A),          \tag{2.1}
\]

then for any fixed

\[
A_2>A_1\ge1/2,
\]

one must have

\[
\boxed{
-\frac1w\log
\frac{P(L\ge\lceil A_2w\rceil)}
     {P(L\ge\lceil A_1w\rceil)}
\to\kappa(p)(A_2-A_1).}                                       \tag{2.2}
\]

The unknown intercept `c_*`, the main-loop cost, and any overall Palm normalization cancel.

In particular the #762 exploratory thresholds `A=1` and `A=2` satisfy the parameter-free exponent prediction

\[
\boxed{
-\frac1w\log
\frac{P(L\ge2w)}{P(L\ge w)}\to\kappa(p).}                     \tag{2.3}
\]

This is a cleaner first test of the mechanism than attempting to identify the complete rate curve from a few widths.

## 3. Finite extra-span ratios inside the linear regime

Suppose a sufficiently uniform local LDP/refined tail asymptotic holds in the saturated regime.  Fix `A>1/2` and an integer `k>=0`.  Since

\[
I_p(A+k/w)-I_p(A)=\kappa k/w,                                  \tag{3.1}
\]

the candidate predicts the local tail ratio

\[
\boxed{
\frac{P(L\ge\lceil Aw\rceil+k)}
     {P(L\ge\lceil Aw\rceil)}
\longrightarrow e^{-\kappa k}}                               \tag{3.2}
\]

up to lattice-periodic/subexponential corrections not determined by the LDP alone.

Thus after the macroscopic bulge has saturated, each additional **microscopic** unit of extreme span costs the ordinary axial mass `kappa`.  A persistent effective rate below `kappa` would exhibit a cheaper network mechanism than the proposed one.

Equation (3.2) is deliberately labelled a refined-tail prediction; (2.2) follows already at the LDP level.

## 4. Morphology implication, stated only at the variational coordinate level

For `A>r_*`, the minimizer of (1.1) remains exactly `r_*` as `A` increases.  Therefore any operational morphology statistic that is rigorously shown to converge to the variational loop coordinate `r` must **saturate** for all `A>=1/2`.

This gives the correct interpretation of the #762 `A=1,2` comparison:

- if a certified winding-core observable is proved to represent the loop coordinate, its `w`-scaled value should be the same at `A=1` and `A=2` to leading order;
- the additional unit of macroscopic span must be carried by the branch part of the optimizer;
- if the same certified core coordinate instead grows proportionally with `A`, the candidate rate mechanism is falsified.

The existing operational `winding core` in #762 is **not yet** proved to be the variational/OZ core.  This note therefore does not equate `L_core` with `2r_*w` or any other particular formula.  The earlier articulation counterexample is precisely why that map must be proved before using morphology to validate the LDP.

## 5. The transition point is itself a directional-norm observable

At differentiability,

\[
2\,\partial_y\tau_p(1,2r_*)=\kappa.                           \tag{5.1}
\]

Thus `r_*` can be predicted from a directional mass certificate without rare-span sampling.  Conversely, once a morphology variable is rigorously tied to `r_*`, it becomes an independent probe of the Wulff shape.

The curvature relation

\[
D^{-1}=\partial_{yy}\tau_p(1,0)                               \tag{5.2}
\]

controls only the small-`A` quadratic start of the rate; `r_*` probes the nonlinear directional norm farther away from the axis.  These are complementary geometric quantities.

## 6. Strong falsification outcomes

The candidate (1.1) should be rejected or revised if any rigorous/asymptotic result shows one of:

1. a large-`A` slope strictly below `kappa`;
2. a nonlinear rate persisting for some `A>=1/2`;
3. a certified loop-coordinate optimizer continuing to grow with `A>=1/2`;
4. a cheaper network topology with the same nonzero deck displacement and span.

A finite-width deviation from (2.3) without controlled errors is not such a falsification.  The useful experiment is an exponent/certificate comparison, not another unconstrained fit.
