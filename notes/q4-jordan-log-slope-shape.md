# Module-specific Q4 Jordan log-slope shape

PR #217 proves that the percolation energy/2-hull Jordan pair survives under
the repository's non-null thermal `Q4` descendant.  Its affine scale cocycle
is necessary but not module-specific: any rank-2 pair at the same exponent
has the same `log L` law.  The first sharper prediction concerns the
coefficient of that logarithm.

For

\[
 Dq=xq,\qquad D\widetilde q=x\widetilde q+q,
\]

a finite dilation gives

\[
 \widetilde q\mapsto s^{-x}(\widetilde q-\log(s)q).
\]

Consequently a top-field response satisfies

\[
 L^{x-2}\delta O_{\rm top}(L,\tau)
 =\widetilde A(\tau)-\lambda_{\rm top}A_q(\tau)\log L,
\]

and, since `log L=(1/2)log N`,

\[
 \boxed{B_{\log N}(\tau)=-\frac{\lambda_{\rm top}}2A_q(\tau)}.
\]

The ordinary torus Ward identity already fixes

\[
 \frac{A_q(\tau)}{A_\epsilon(\tau)}=\frac{493}{96}g_2(\tau),
\]

so before the nonuniversal common overlap the coefficient is
`-493/192`.  If the residual-to-slope/root normalization removes the same
thermal-primary block as in the existing closure, then

\[
 B_{\rm root}(\tau)=C_J\operatorname{Re}g_2(\tau)
\]

(or `C_J g2` for a chiral projection).  This freezes

\[
 B_{\rm root}(e^{i\pi/3})=0,\qquad
 \frac{B_{\rm root}(2i)}{B_{\rm root}(i)}=\frac{11}{4}
\]

in the area-normalized convention.  The companion script independently
checks the CM zero and rectangular ratio by a 90-decimal direct E4 q-series,
in addition to consuming the repository's exact Ward and Hecke fractions.

The score must be covariance aware: log slopes, bottom-Q4 amplitudes, root
normalization, and final ratios or residuals are recomputed inside each common
replicate.  A useful amplitude-free cross-channel quantity is

```text
[B_Sprime(tau1)/A_bottom(tau1)] / [B_Sprime(tau2)/A_bottom(tau2)] = 1.
```

## Boundary

This is not a derivation of the complete top torus one-point.  The additive
shape `A_tilde(tau)` changes under `q_tilde -> q_tilde+alpha q` and also
depends on generic-c derivatives of the colliding energy and hull blocks.
No unique `E2`, `log eta`, weight derivative, or modular derivative is fixed
by the Jordan relation alone.  Likewise, a norm-4 affine scale cocycle alone
cannot identify the Q4 Jordan module; module specificity requires the joint
logarithmic scale law and the independently reconstructed bottom-Q4 shape.

Reproduce with:

```bash
python scripts/derive_q4_jordan_log_slope_shape.py --dps 90
python -m unittest discover -s tests -p 'test_q4_jordan_log_slope_shape.py'
```
