# P398: the existing pair density sees a strongly ray-selective current

**Yes.** The single predefined pair consisting of the source and size-two
cluster charge detects the orientation invisible to real self-correlations.
It captures 53.1666% of the minus current-source squared norm but only
0.433293% of the plus norm. This is not a small total current in plus:
the full normalized current-source norms squared are almost equal,
0.3808266892 and 0.3600854185.

Parent: `520a9d21`. Definition:
`analysis/p398_width8_pair_density_current.json`. Executable and complete
fixed-grid score: `scripts/p398_width8_pair_density_current.py` and
`results/p398-width8-pair-density-current/latest.json`.

## Fixed pair and phase convention

Use the same site origin, character i, and protected rays as before:

\[
\psi_s=(A+s\zeta L)/\sqrt2,\quad \zeta=e^{-i\pi/4},\qquad
t_s=P_s T_2,\quad T_2=\sum_{|C|=2}\sum_{j\in C}i^j.
\]

P_s is the orthogonal stationary-L2 projector onto the same Kreweras ray.
Let e=psi_s/||psi_s|| and

\[
\eta_s=\frac{t_s-e\langle e,t_s\rangle_\pi}
 {\|t_s-e\langle e,t_s\rangle_\pi\|_\pi}.
\]

The denominator is positive real: **there is no extra rephasing**. Both the
raw normalized pair (e,t_s/||t_s||) and orthonormal pair (e,eta_s) are saved.
Every function has zero stationary mean because its cyclic character is i.

For C_fg(t)=<f,exp(tG)g>_pi, report the anti-Hermitian **difference**, not
half the difference:

\[
\mathcal A_{fg}(t)=C_{fg}(t)-\overline{C_{gf}(t)},\qquad
\mathcal A'_{fg}(0)=2\langle f,Jg\rangle_\pi.
\]

For G* the entire quantity changes sign by the adjoint identity. For S it
is identically zero. Neither control was solved again.

## Initial direction and captured current

| Quantity | Minus ray | Plus ray |
|---|---:|---:|
| Var(t_s) | .4361178135 | .4514707838 |
| Var(t_s orthogonal to source) | .2119873684 | .3883270749 |
| ||J e||^2 | .3808266892 | .3600854185 |
| |<eta,J e>|^2 | .2024725145 | .0015602245 |
| Fraction of ||J e||^2 captured | **.5316657688** | **.0043329289** |
| Raw normalized anti-cross derivative, coefficient of 1+i | -.4436608580 | +.0518074909 |
| Orthonormal anti-cross derivative, coefficient of 1+i | **-.6363529122** | **+.0558609798** |

The initial directions are opposite in this fixed physical gauge, and the
orthonormal derivative magnitude is about 11.39 times greater on minus.
The phase of an individual cross-observable changes if its readout is
rephased; the capture fraction, nonzero status, and a sign reversal versus
distance are invariant under any constant partner rephasing. These are
projections of a current, not independent definitions of two global current
orientations.

There is a useful exact connection to the existing memory calculation.
The established geometry identity gives

\[
G\psi_s=-3\psi_s+\sqrt2\,s\zeta t_s.
\]

Thus eta_s is exactly the previously named first forward-force direction,
up to its explicitly fixed phase. Write M=-G, c=(I-ee*)Me, F=||c||^2, and
k_0=<e,M(I-ee*)Me>. Since the present self-correlations are real,

\[
\langle c,Je\rangle=(k_0-F)/2,\qquad
\boxed{\frac{|\langle\eta_s,Je\rangle|^2}{\|Je\|^2}
 =\frac{(k_0-F)^2}{4F\|Je\|^2}.}
\]

For minus, k_0/F=3.50896; for plus, k_0/F=1.04607. The plus pair-density
force has nearly the same forward variance and left-right feedback, so
their current-sensitive difference is small, even though both feedback
and the full current norm are substantial. The earlier large plus memory
does **not** imply that this particular pair reads a large oriented current.

## A lag-dependent sign reversal appears without searching a new lag

All entries below give A_e,eta(s) divided by (1+i), in the predefined grid.

| s | Minus | Plus |
|---|---:|---:|
| 0 | 0 | 0 |
| .05 | -.0267960327 | +.0017821829 |
| .1 | -.0451303768 | +.0021309949 |
| .25 | -.0673956881 | -.0000538937 |
| .5 | -.0573160321 | -.0033162175 |
| 1 | -.0211961465 | -.0028350833 |
| 2 | -.0015985015 | -.0004892052 |
| 4 | -.0000060994 | -.0000099007 |

Minus is negative at every sampled positive distance. Plus is positive at
.05 and .1, then negative from .25 onward: an ordered pair-correlation
reversal is already bracketed by the fixed grid. No extra root search was
performed. This is a lag-dependent projection of a *fixed stationary*
generator, not a change in the underlying stationary probability current.
At s=4 the absolute plus projection even exceeds minus despite the much
smaller initial derivative. A single instantaneous pair-current amplitude
cannot encode the full distance-dependent directional propagation.

## Scope

Only {psi_s,P_s T2} was used; no alternative readout, width, fit, simulation
or test campaign was added. This deterministic result shares the exact same
finite-process dependency block as the memory and reversible-control
results. It identifies an observable of irreversibility inside that process,
not a site-Matching field, continuum spin, Jordan mode, or intrinsic
morphism-history effect.
