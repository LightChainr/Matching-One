# P398: T4 repairs hidden transport, not a missing initial current direction

The T2 pair's plus blindness was an **observer-direction issue**, not a
weak current. But correcting that observer does not explain T4's tail repair.
The existing 7-dimensional plus geometry already captures 81.7318% of the
current-source squared norm; T4 adds only .16723 percentage points.
Keeping the entire initial current direction in a two-observable model
still loses nearly all of the plus long tail.

This is one generator-determined readout, not another geometric search.
The definition and executable are
`analysis/p398_width8_current_source_geometry.json` and
`scripts/p398_width8_current_source_geometry.py`; all scores are in
`results/p398-width8-current-source-geometry/latest.json`.

## The unique direction and what is tautological

With the unchanged normalized source e=psi_s/||psi_s||, define

\[
r_s=\frac{Je-e\langle e,Je\rangle_\pi}
 {\|Je-e\langle e,Je\rangle_\pi\|_\pi},\qquad J=(G-G^*)/2.
\]

The denominator is positive real; the psi phase is unchanged and r is not
rephased. Here <e,Je>=0, so r=Je/||Je||. The initial anti-Hermitian cross
derivative is therefore -2||Je||: -1.2342231389 on minus and -1.2001423557
on plus. This saturates the Cauchy bound **by construction**; it is not a
discovery of a specially fitted optimal observable.

It nevertheless establishes that the .433% T2 capture on plus was a
direction mismatch. No geometric interpretation is automatically assigned
to r: it uses the known stationary time reversal, hence pi-weighted
transition structure.

## How much do the existing geometries contain?

| Fraction of ||Je||^2 | Minus | Plus |
|---|---:|---:|
| 7-dimensional span | 93.22950% | 81.73183% |
| 8-dimensional span | 93.35889% | 81.89906% |
| Added by T4 | **.12939 percentage points** | **.16723 percentage points** |

The increment is computed both from the projector difference and from the
overlap with the single normalized T4 residual orthogonal to the old
7-space. No new rank or columns are selected.

This tiny direct gain accompanies the previously measured 5.6856-point
plus t=4 repair in G, which disappears in S. T4 is thus not principally
rescuing the *initial* current-source direction. Its important effect is
consistent with an indirect bridge through hidden geometric propagation.
The projection numbers do not, by themselves, identify every edge of that
bridge or establish a unique microscopic pathway.

## Capturing the whole current still fails as a propagation closure

Keep the same orthonormal pair (e,r) for both G and S. The original mass
matrices M=-G are approximately

\[
B_- =\begin{pmatrix}3.36882024&1.49576347\\
 .26154033&3.77897026\end{pmatrix},\qquad
B_+ =\begin{pmatrix}3.69141527&1.31300846\\
 .11286611&5.98235033\end{pmatrix}.
\]

For S simply take the Hermitian part, keeping r fixed from the original G.

| Model | Minus slow mass | Plus slow mass | Minus t=4 error | Plus t=4 error |
|---|---:|---:|---:|---:|
| Full G reference | 2.81965863 | 1.95575014 | — | — |
| (e,r) under G | 2.91567171 | 3.62845820 | -23.6382% | **-99.7429%** |
| Full S reference | 2.54079598 | 1.83631805 | — | — |
| Same (e,r) under S | 2.67162871 | 3.48766990 | -31.5462% | **-99.7268%** |

Neither two-dimensional projection has a ray crossing in the original
.05..1 bracket. These are observable projections, not new Markov state
chains. Complete correlations at every original fixed lag are saved.

There is an exact explanation for the lost initial feedback, not just a
bad fit. Let P project onto (e,r), Q=I-P, and k0 denote initial logarithmic
curvature. Since QJe=0,

\[
\boxed{k_{0,G}^{\rm full}-k_{0,G}^{P}
=\langle e,GQGe\rangle
=\|QSe\|^2
=k_{0,S}^{\rm full}-k_{0,S}^{P}.}
\]

Indeed QGe=QG*e=QSe. The missing curvature is exactly **.0602551994** on
minus but **2.9272489452** on plus. The plus two-dimensional G projection
retains only .1481941519 of the full 3.0754430971 initial curvature; its S
counterpart retains .5082795704 of 3.4355285156.

After retaining the complete instantaneous current direction, the omitted
initial feedback is therefore *strictly reversible-force geometry*. A
current-sensitive observation alone is not a sufficient dynamical state.
Combined with the preceding T4 counterfactual, the result separates the
initial current readout from the subsequent current-dependent mixing of
geometric memory.

## Full ordered cross-correlation, without projecting away hidden states

For the same (e,r), the full anti-cross C_er-conj(C_re) is real in this
fixed phase convention:

| Distance | Minus | Plus |
|---|---:|---:|
| .05 | -.05172559 | -.04798995 |
| .1 | -.08675495 | -.07716526 |
| .25 | -.12828842 | -.10332006 |
| .5 | -.10795711 | -.07994826 |
| 1 | -.03950763 | -.03127015 |
| 2 | -.00295886 | -.00430601 |
| 4 | -.00001126 | -.00008576 |

Both rays are negative throughout the fixed grid, unlike the plus T2 pair's
lag reversal. The late plus excess is visible in the *full* cross-response
but absent from the two-dimensional closure. No additional zero search,
lag optimization, simulation, width expansion, or repeated test campaign
was performed. All results remain one deterministic finite-process block,
not an independent field/Jordan/morphism-history identification.
