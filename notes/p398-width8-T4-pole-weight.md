# P398: T4 repairs the plus tail by accelerating its slow pole

**The principal repair is pole movement, not removal of spectral weight.**
The plus slow mass increases from 1.9313697560 to 1.9479283947; its source
residue also *increases*, partially opposing the reduction in the t=4 tail.
The signed t=4 correction is -5.68562 percentage points of the full-G
reference: -6.96791 from the pole, +1.28154 from weight, and +.000753 from
all remaining modes.

Only the saved seven/eight-dimensional original-G Schur blocks from
`074a5f53` were used. There was no reconstruction or solve of the 1430-state
process, new lag, simulation, or coupling fit. The script is
`scripts/p398_width8_T4_pole_weight.py`; the complete score is
`results/p398-width8-T4-pole-weight/latest.json`.

## Separate the two effects without refitting

Write the slow contribution as L(t)=r exp(z t), where z=-mu is the
generator pole and r its normalized source residue. Both endpoint residues
are positive. The exact multiplicative separation is

\[
\log(L_8/L_7)=t(z_8-z_7)+\log(r_8/r_7).
\]

For an additive comparison, use the explicitly symmetric two-factor split

\[
\Delta L=\underbrace{\frac{r_7+r_8}{2}(e^{z_8t}-e^{z_7t})}_{\text{pole}}
+\underbrace{\frac{e^{z_7t}+e^{z_8t}}2(r_8-r_7)}_{\text{weight}}.
\]

The remaining change in the total correlation is retained as the other-mode
term. This is an exact output decomposition under the stated convention,
not an additive microscopic causal percentage.

| Ray / model | Slow mass mu | Slow source residue r |
|---|---:|---:|
| Minus, old 7 | 2.8428902867 | .6161898229 |
| Minus, new 8 | 2.8404109327 | .6140295644 |
| Minus, full G | 2.8196586326 | .5838283762 |
| Plus, old 7 | 1.9313697560 | .4612922882 |
| Plus, new 8 | 1.9479283947 | .4669440563 |
| Plus, full G | 1.9557501384 | .4706022096 |

Both mass and residue move toward their respective full-G values on each
ray. But their effects on the scalar tail oppose each other: correcting a
residue parameter need not itself reduce the current correlation error.

At the original t=4, in percentage points of that ray's full-G correlation:

| Ray | Pole part | Weight part | Other modes | Total change |
|---|---:|---:|---:|---:|
| Minus | +.952540 | -.337320 | +.004680 | **+.619901** |
| Plus | -6.967914 | +1.281542 | +.000753 | **-5.685620** |

Under this additive convention, pole movement supplies 122.55% of the net
plus repair and weight change offsets 22.54%; the other-mode term offsets
only .01325%. On minus, pole slowing supplies 153.66% of the net uplift,
while the residue reduction offsets 54.42%.

The convention-independent slow-contribution log changes agree:

| Ray | Pole log change | Weight log change | Total slow log change |
|---|---:|---:|---:|
| Minus | +.0099174160 | -.0035119926 | +.0064054234 |
| Plus | -.0662345548 | +.0121775838 | -.0540569710 |

The bridge therefore does not mainly siphon plus spectral weight away. It
corrects an artificially slow reduced pole while *restoring* some source
weight that had previously partially cancelled that pole error.

## What the Schur denominator says

Retain the previous block and phase convention H8=[[A,b],[c,d]]. Introduce
only an analytic feedback bookkeeping parameter,

\[
H(\lambda)=\begin{pmatrix}A&b\\\lambda c&d\end{pmatrix},
\qquad D(z,\lambda)=z-d-\lambda c(z-A)^{-1}b.
\]

At lambda=0 its source correlation is the old seven-state projection; at
1 it is the new eight-state projection. Intermediate matrices are **not**
asserted to be positive Markov generators. No intermediate value was
scanned or selected.

For a simple new pole away from the old poles, R=(z-A)^(-1) gives

\[
z'(\lambda)=\frac{cRb}{1+\lambda cR^2b},\qquad
r(\lambda)=\frac{\lambda(e^*Rb)(cRe)}{1+\lambda cR^2b}.
\]

At lambda=1, the denominator derivatives are **+2116.6161** on minus and
**-381.3972** on plus. Although cRb=z-d is positive on both, these opposite
signed feedback slopes move their generator poles in opposite directions.
The residues recovered from the Schur factors are .6140295644 and
.4669440563, matching the source spectral projectors.

This is the source-visible action of the non-selfadjoint bridge, not a new
slow bare T4 mode. The bare d poles remain near -9.3; signed geometric
feedback shifts a much slower existing mode.

Analytic eigenprojector differentiation at only the two endpoints gives:

| Ray / endpoint | d mu / d lambda | d r / d lambda | Pole part of d log L | Weight part of d log L |
|---|---:|---:|---:|---:|
| Minus, 0 | -.0020161397 | -.0016351587 | +.0080645586 | -.0026536606 |
| Minus, 1 | -.0030482265 | -.0028385650 | +.0121929061 | -.0046228475 |
| Plus, 0 | +.0142398366 | +.0044306620 | -.0569593463 | +.0096048906 |
| Plus, 1 | +.0192743653 | +.0071732450 | -.0770974612 | +.0153621081 |

Both endpoints show the same local pole-dominated mechanism. No global
monotonicity or absence-of-collision claim over uncomputed intermediate
lambda values is needed for the endpoint decomposition.

All of this concerns the same fixed original-G projection. The preceding
S result already establishes that the main plus repair is current-enabled;
no new S campaign, low-order null analysis, continuum field claim, or
independent evidence is added here.
