# P398: fast-to-slow inversion survives deletion of every stationary current

The answer is **no: irreversible circulation is not required for this
inversion**. The one predefined counterfactual retains the same 1430 states,
stationary geometry, AP/landing sources, source variances, and individual-state
exit rates. It replaces the backward generator G by S=(G+G*)/2, with no fit,
parameter interpolation, width change, or Monte Carlo.

Definition commit: `5e47bdb6`; parent geometry result: `30eef34a`.
Executable: `scripts/p398_width8_reversible_current_control.py`.
Full fixed-grid outputs, spectra, residues and input hashes:
`results/p398-width8-reversible-current-control/latest.json`.

## The propagation result

The original Kreweras map preserves pi and commutes with G. It is consequently
unitary in stationary L2, also commutes with G*, and protects the same two
93-dimensional rays under S. Nothing is reidentified between the processes.

| Quantity | Original G | Current-deleted S |
|---|---:|---:|
| Initial decay, minus | 3.3688202414 | 3.3688202414 |
| Initial decay, plus | 3.6914152676 | 3.6914152676 |
| Lowest visible mass, minus | 2.8196586326 | 2.5407959794 |
| Lowest visible mass, plus | 1.9557501384 | 1.8363180504 |
| Normalized ray crossing | 0.2656573200 | 0.2722634760 |
| Initial log curvature, minus | 0.4514576654 | 0.8322843545 |
| Initial log curvature, plus | 3.0754430971 | 3.4355285156 |
| Integrated projected memory, minus | 0.1212071174 | 0.2291667011 |
| Integrated projected memory, plus | 0.7859741072 | 0.8811456488 |

Both processes have a plus source with **larger initial decay and smaller
lowest visible mass**. The crossing shifts by only +2.4867%. Under S the
leading normalized residues are 0.53072405 (minus) and 0.45784111 (plus),
and all spectral weights are nonnegative: the reversal already occurs in
ordinary reversible positive mixtures of exponential modes. It does not
need negative residues, complex modes, or Jordan structure.

The normalized plus/minus ratio on the fixed lag grid is:

| Distance s | Original G | Current-deleted S |
|---|---:|---:|
| .05 | .98719822 | .98716658 |
| .1 | .98074746 | .98059499 |
| .25 | .99570258 | .99415226 |
| .5 | 1.11946659 | 1.10803151 |
| 1 | 1.68824298 | 1.59004251 |
| 2 | 4.31451524 | 3.43824768 |
| 4 | 25.41893805 | 14.42424134 |

All seven individual correlations increase after current deletion in this
example. That observation is not asserted as a universal pointwise ordering.
The circulation nevertheless increases the **relative** late separation:
the lowest-mass gap rises from 0.7044779290 to 0.8639084942, a 22.631% increase.
Thus the current is an amplifier of the two-ray contrast, not its source.

## What the current does at short distance

The prewritten identity gives, on either ray,

\[
\kappa_S(0)-\kappa_G(0)
=\frac{\|J\psi\|_\pi^2}{\|\psi\|_\pi^2},\qquad
\kappa=(\log u)''.
\]

The right sides are **0.3808266892** for minus and **0.3600854185** for plus.
They match the curvature changes to below 6e-15 in the same matrix readout.
In reversible dynamics the curvature is the variance of the source-visible
mass distribution. The skew part suppresses that curvature by 45.757% on
minus but only 10.481% on plus. It raises the plus/minus initial-memory ratio
from 4.12783 to the previously reported 6.81225.

Crucially, the two absolute current-force contributions are nearly equal.
Deleting them changes the **difference** of initial curvatures by only
-0.79045%, while leaving the difference of initial slopes exactly fixed.
This explains why the early crossing barely moves even though each ray's
memory and long tail change appreciably. The unchanged weak-source
variance cancellation and reversible source/mode overlaps are already
sufficient for the inversion; nonreversibility adds differential long-time
acceleration, stronger on the minus ray.

There is also an exact limit to what self-correlations can identify:

\[
\operatorname{Re}\langle f,e^{tG}f\rangle_\pi
=\operatorname{Re}\langle f,e^{tG^*}f\rangle_\pi.
\]

Reversing every current therefore leaves the real self-correlation unchanged.
The present G-to-S intervention identifies an effect of *current magnitude*,
not its orientation. For a named multi-observable correlation matrix C(t),
the directional object is its anti-Hermitian part:

\[
\left.\frac{d}{dt}[C_{fg}(t)-\overline{C_{gf}(t)}]\right|_{t=0}
=2\langle f,Jg\rangle_\pi.
\]

Such a test must use observables for which symmetry does not force this
cross term to vanish, for example a source and a suitable named geometric
partner within the same protected ray. No new cross-correlation selection
or calculation is included here.

## Meaning and boundary

This is a deterministic finite-process intervention, not new statistical
evidence from independent samples. S is a genuine reversible generator:
nonnegative off-diagonal rates, stationary residual below 9e-16,
detailed-balance residual below 4e-18, and exactly unchanged diagonal in the
numeric representation. No broad regression suite was rerun; these are the
arithmetic identities of this one intervention.

The pi-weighted reverse edges need not be the original local square-bond
row word. The result establishes sufficiency of reversible dynamics on the
same stationary connectivity geometry, **not** a second microscopic lattice
model, a continuum field count, or a site-Matching/Jordan identification.
The named geometric compression and current-deletion answer different
questions on the same underlying finite process and must not be counted as
independent evidence.
