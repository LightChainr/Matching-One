# P398: the defined stationary-current deletion

This proposition and comparison are written before the new numerical output.
Let G be the backward row generator, with nonnegative off-diagonal rates,
G1=0, and a strictly positive stationary column pi satisfying G^T pi=0.
Write Pi=diag(pi), and define

\[
G^*=\Pi^{-1}G^T\Pi,\qquad S=(G+G^*)/2,\quad J=(G-G^*)/2.
\]

Then S is a valid reversible generator: its off-diagonal rates are averages
of forward and reversed rates, S1=0, and Pi S is symmetric. It preserves pi
and all stationary readout variances. Since G* has the same diagonal as G,
S also preserves the total exit rate at every state. Strict positivity of pi
is needed for this formula; otherwise first restrict to its invariant support.
J is skew-adjoint in stationary L2 and generally is not a Markov generator.

For a fixed source f, the real Dirichlet initial slope is unchanged. If
`<f,Jf>_pi=0`, as for the real-valued ray correlations in this comparison,

\[
\boxed{u_S''(0)-\operatorname{Re}u_G''(0)
=\|Jf\|_\pi^2/\|f\|_\pi^2.}
\]

The same formula holds for the real logarithmic curvature when the initial
slopes are real. For a general complex source, the logarithmic version
subtracts `|<f,Jf>_pi|^2/||f||_pi^4` from the right side. Thus deleting
currents does not necessarily remove memory: it increases this initial
curvature. A reversible projected memory kernel is a positive mixture
`c^* exp(-D t)c` in the orthogonal stationary projection.

The only comparison is G versus S on the same fixed psi-minus/plus, same
lag grid .05,.1,.25,.5,1,2,4, crossing bracket .05..1, and lowest visible
masses. No current-strength interpolation or source refit will be tried.

S preserves the stationary geometry but can introduce pi-weighted reverse
transitions not expressible as the original local bond update word. It is a
controlled finite Markov process, not automatically a new square-bond lattice
transfer model or an equilibrium continuum field identification.
