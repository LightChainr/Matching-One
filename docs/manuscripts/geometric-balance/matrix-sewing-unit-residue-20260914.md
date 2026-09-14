# Matrix sewing lemma: finite memory preserves the unit logarithmic residue

2026-09-14.  This is a self-contained analytic lemma for a finite-state Markov-additive/cyclic kernel.  It generalizes the scalar renewal-loop calculation already on the branch and isolates what **cannot** be blamed on finite local memory in the still-missing SITE sewing theorem.

## 1. Finite-state cyclic kernel

Let

\[
A(z,y)=\sum_{x\ge1}\sum_{j\in\mathbb Z}A_{x,j}z^xy^j
\]

be a `d x d` matrix Laurent polynomial with nonnegative real matrices `A_{x,j}` and finite support.  The longitudinal increment `x` is strictly positive.  Assume:

1. for real `z>0,y=1`, the nonnegative matrix `A(z,1)` is irreducible in a neighbourhood of the critical point;
2. there is a unique `R>0` with Perron root
   \[
   \rho(A(R,1))=1;
   \]
3. this Perron eigenvalue is simple, all other eigenvalues of `A(R,1)` have modulus strictly below one, and the joint support is aperiodic so `(R,1)` is the only dominant singular point modulo the Fourier period;
4. transverse reflection symmetry gives zero mean transverse increment under the critical Perron tilt and an effective variance `sigma_eff^2>0`.

Put

\[
\kappa=\log R.
\]

Define the cyclic loop coefficient

\[
\boxed{
L_w=w[z^wy^0]\{-\log\det(I-A(z,y))\}.}                         \tag{1.1}
\]

The logarithm has the exact formal expansion

\[
-\log\det(I-A)=\sum_{n\ge1}\frac1n\operatorname{tr}A^n.       \tag{1.2}
\]

Thus (1.1) is the finite-memory analogue of the scalar cyclic-renewal object: the trace closes the internal state, `1/n` unmarks the cyclic transition index, and `w` restores longitudinal translation.

## 2. Perron transform and the diffusion constant

Let `r,l` be positive right/left Perron vectors of `A(R,1)`, normalized by `l^T r=1`.  The critical tilted transition on an edge carrying increment `(x,j)` is proportional to

\[
(A_{x,j})_{ab}R^x\frac{r_b}{r_a}.
\]

After normalization this is a finite-state Markov-additive chain.  Let

\[
\mu=E X>0
\]

be its mean longitudinal advance and let `sigma_eff^2` be the asymptotic variance per transition of the accumulated transverse displacement, including the state correlations.  Set

\[
\boxed{D=\sigma_{eff}^2/\mu.}                                  \tag{2.1}
\]

Equivalently, let `lambda(z,t)` be the analytic Perron eigenvalue of `A(z,e^t)` near `(R,0)` and write `z=Re^u`.  Standard analytic perturbation of a simple eigenvalue gives

\[
\partial_u\log\lambda(0,0)=\mu,
\qquad
\partial_t\log\lambda(0,0)=0,
\qquad
\partial_t^2\log\lambda(0,0)=\sigma_{eff}^2.                  \tag{2.2}
\]

The last derivative is the Green--Kubo/asymptotic variance of the additive functional, not merely a one-step variance when internal states carry memory.

## 3. The unit-residue theorem

**Theorem.**  Under the assumptions above,

\[
\boxed{
L_w=\frac{e^{-\kappa w}}{\sqrt{2\pi D w}}
\left(1+O(w^{-1})\right).}                                    \tag{3.1}
\]

In particular, a finite amount of Markov memory changes `kappa` and `D` but creates **no independent multiplicative residue** in the normalized cyclic object (1.1).

### Proof

Factor the determinant locally using the simple Perron band:

\[
\det(I-A(z,e^{i\theta}))
=(1-\lambda(z,i\theta))G(z,\theta),                            \tag{3.2}
\]

where `G` is analytic and nonzero near `(R,0)`.  Hence

\[
-\log\det(I-A)
=-\log(1-\lambda)+\text{analytic}.                             \tag{3.3}
\]

The analytic term has a strictly larger radius in the dominant direction (after shrinking the neighbourhood and using the assumed global aperiodicity), so it is exponentially/subdominantly smaller than the Perron singularity in the coefficient under discussion.

For small real `theta`, the implicit function theorem supplies the unique root `R(theta)` near `R` of

\[
\lambda(R(\theta),i\theta)=1.                                  \tag{3.4}
\]

Write `R(theta)=R exp u(theta)`.  Expanding `log lambda` with (2.2) at imaginary transverse tilt gives

\[
0=\mu u(\theta)-\frac12\sigma_{eff}^2\theta^2+O(\theta^4),
\]

therefore

\[
\boxed{
\log R(\theta)=\log R+\frac12D\theta^2+O(\theta^4).}           \tag{3.5}
\]

Near its simple zero,

\[
1-\lambda(z,i\theta)
=C(\theta)(1-z/R(\theta))+O((1-z/R(\theta))^2)
\]

with `C(theta) !=0`.  The constant `C(theta)` enters only the **analytic** additive term `-log C(theta)`.  The logarithmic singular part is exactly

\[
-\log(1-z/R(\theta)),
\]

whose coefficient satisfies the identity

\[
w[z^w]\{-\log(1-z/R(\theta))\}=R(\theta)^{-w}.                 \tag{3.6}
\]

This is the source of the unit residue: Perron left/right overlaps and the derivative of the eigenvalue do not multiply the logarithmic coefficient.

Fourier inversion in the transverse displacement now gives

\[
L_w=\frac1{2\pi}\int_{-\pi}^{\pi}R(\theta)^{-w}\,d\theta
+\text{subdominant}.                                           \tag{3.7}
\]

Aperiodicity makes the integral away from `theta=0` exponentially smaller.  Insert (3.5) in a neighbourhood of zero and apply the ordinary one-dimensional Laplace method:

\[
\frac{R^{-w}}{2\pi}\int_{\mathbb R}
 e^{-Dw\theta^2/2}\,d\theta
=\frac{R^{-w}}{\sqrt{2\pi Dw}}.
\]

The finite-support analytic expansion gives the stated `O(1/w)` correction.  Since `R^{-w}=e^{-kappa w}`, (3.1) follows.  `square`

## 4. Directional curvature is automatically the inverse diffusion constant

The same Perron cumulant also controls the open-endpoint large-deviation cost of the Markov-additive chain.  For endpoint transverse slope `v`, the local rate expansion is the Legendre transform of the transverse cumulant per unit longitudinal distance:

\[
\tau(1,v)=\kappa+\frac{v^2}{2D}+O(v^4).                        \tag{4.1}
\]

Therefore

\[
\boxed{\partial_{vv}\tau(1,0)=D^{-1}.}                        \tag{4.2}
\]

Writing the homogeneous norm as

\[
\tau(r\cos\theta,r\sin\theta)=r\kappa(\theta)
\]

and using reflection symmetry gives

\[
\boxed{D^{-1}=\kappa(0)+\kappa''(0).}                          \tag{4.3}
\]

Thus the diffusion/Wulff-curvature relation highlighted in `structural-consequences-20260914.md` is not a separate miracle in a finite-memory sewing model: it is an automatic consequence of the same twisted Perron band.

## 5. What this settles about `sewing-with-memory.md`

The actual complete SITE activity has a three-column local weight and therefore microscopic memory.  This lemma shows:

> **Finite local memory, by itself, does not change the `w^{-1/2}` closure power and does not create an arbitrary residue.**

If an exact percolation sewing theorem identifies complete components with (1.1) for a finite-state kernel satisfying the hypotheses, then automatically

\[
\beta=1/2,\qquad \zeta=1,
\]

in the notation of the branch's HK/unit-residue conjectures.

Therefore any non-unit `zeta` or different power for the actual component density must come from a failure of at least one of the following identifications:

1. the complete-component class is not exactly the cyclic trace/log-determinant object;
2. the natural mark/cut law introduces an extra insertion rather than pure cyclic unrooting;
3. the effective state is genuinely infinite and lacks a simple isolated Perron band/quasi-compact reduction;
4. more than one soft transverse band contributes at the same exponential rate;
5. lattice periodicity/aperiodicity leaves multiple dominant saddles.

This is much narrower than saying that the three-column boundary memory itself leaves the prefactor unknown.

## 6. Countable-state extension: the actual next theorem

The small-`p` programme in `research-frontier-20260914.md` should aim to replace the finite matrix by a positive operator on a weighted Banach space of column/connectivity states.  A sufficient package would be:

- quasi-compactness and an isolated simple Perron eigenvalue near the physical root;
- analytic dependence on longitudinal/transverse fugacities;
- an aperiodic Markov-additive Perron transform with finite exponential moments;
- trace/Fredholm-determinant control strong enough that the cyclic connected object has the same single logarithmic singularity.

Under those hypotheses the proof above carries over with `det` replaced by a suitable Fredholm determinant.  Establishing those operator hypotheses for actual SITE complete components on a nonempty small-`p` interval would close the prefactor problem there without importing a two-point OZ amplitude.

## 7. Claim boundary

The finite-matrix theorem is an analytic statement about the explicitly defined cyclic object (1.1), not yet an identification of that object with actual SITE complete components.  It does, however, remove **finite Markov memory** from the list of possible reasons for an anomalous power/residue and supplies the exact curvature--diffusion relation once the common kernel is established.
