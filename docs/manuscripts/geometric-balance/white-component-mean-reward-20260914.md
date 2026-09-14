# Mean volume and boundary rewards of the giant complementary white component

2026-09-14.  This note upgrades part of `supercritical-white-slab-bulk-20260914.md` from conjecture to an author-level theorem.  The **samplewise slab LLN** remains open, but the component-Palm means are fixed by a local essential-membership limit plus exact mass transport.

## 1. Setup

Fix black NN site density

\[
0<p<p_c(G4),
\qquad
q=1-p>p_c(G8).                                                \tag{1.1}
\]

Work on the cylinder

\[
C_w\times\mathbb Z.                                           \tag{1.2}
\]

Black uses NN connectivity and white uses matching connectivity from the same Bernoulli labels.

Let

\[
\nu_w=\nu_w^4(p)=\nu_w^8(q)                                  \tag{1.3}
\]

be the common row intensity of complete essential components.  For a white essential component under component Palm, write

\[
L_w=\text{vertical span},
\quad
N_w=\text{number of white sites},
\quad
B_w=\text{number of distinct black external boundary sites}. \tag{1.4}
\]

The previous alternating-barrier theorem gives

\[
\nu_wE_WL_w\to1.                                             \tag{1.5}
\]

## 2. Probability that a fixed cylinder site belongs to a white essential component

Let

\[
\rho_w(q)
=P_q^{C_w\times\mathbb Z}
(0\text{ belongs to a white matching essential component}).   \tag{2.1}
\]

Let

\[
\theta_8(q)
=P_q^{\mathbb Z^2}(0\leftrightarrow\infty\text{ in }G8).      \tag{2.2}
\]

We claim

\[
\boxed{\rho_w(q)\to\theta_8(q).}                              \tag{2.3}
\]

### 2.1 Upper bound

If the cylinder component of the origin is horizontally essential, lift it to the plane.  The lift contains a path from the origin to a nonzero horizontal translate of itself, hence it exits every Euclidean ball of radius `<w/2`.

Take `R=floor(w/4)`.  The ball `B_R` injects into the cylinder, so its SITE variables have exactly the planar iid law.  Therefore

\[
\rho_w(q)
\le P_q^{\mathbb Z^2}(0\leftrightarrow\partial B_R).          \tag{2.4}
\]

As `R->infinity`, the decreasing one-arm events converge to membership in the infinite cluster.  Hence

\[
\limsup_{w\to\infty}\rho_w(q)\le\theta_8(q).                 \tag{2.5}
\]

### 2.2 Lower bound by a high-probability white frame

Because black NN density `p` is subcritical, black crossings of a rectangle in its long direction have probability at most `C e^{-cw}` uniformly over any fixed finite collection of aspect ratios.

For the square-site matching pair, the standard finite-rectangle matching dichotomy says that failure of a white matching crossing in one direction forces a black NN crossing in the transverse direction.  Consequently all white matching rectangle crossings in a fixed finite gluing scheme of scale `w` occur with probability

\[
1-O(e^{-cw}).                                                  \tag{2.6}
\]

Choose constants

\[
0<r_1<r_2<1/4.                                                \tag{2.7}
\]

Using finitely many such high-probability crossings, construct a connected white **frame** `F_w` with the following deterministic properties:

1. inside the injected annulus `B_{r_2w}\setminus B_{r_1w}`, it contains a closed white matching circuit surrounding the origin;
2. this circuit is connected, through a fixed number of white crossing rectangles, to a white horizontal winding ring on the cylinder;
3. the whole construction uses a vertical band of height `O(w)`.

The gluing is the supercritical version of the staircase/corridor constructions already used in the geometric manuscript.  Since the number of rectangles is fixed,

\[
P(F_w^c)\le C e^{-cw}.                                        \tag{2.8}
\]

On `F_w`, any white path from the origin to `partial B_{r_2w}` must intersect the surrounding white circuit in the embedded matching spine, and hence joins the horizontal winding ring.  Therefore its cylinder component is essential.

The ball `B_{r_2w}` injects, so

\[
\rho_w(q)
\ge
P_q^{\mathbb Z^2}(0\leftrightarrow\partial B_{r_2w})
-Ce^{-cw}.                                                     \tag{2.9}
\]

Letting `w->infinity`,

\[
\liminf\rho_w(q)\ge\theta_8(q).                              \tag{2.10}
\]

Together with (2.5), this proves (2.3).

The only planar input beyond subcritical black sharpness is the finite matching crossing dichotomy already inherent in the square-site 4/8 convention; no supercritical white correlation length is imported.

## 3. Exact Campbell identity for white component volume

Anchor every white essential component at its lowest row as in the preceding notes.  For each component `C` with anchor row `a(C)`, send one unit of mass from `a(C)` to every white site of `C`, recording the target row.

Stationarity in the vertical coordinate gives the discrete Campbell/mass-transport identity

\[
\boxed{
\nu_w E_WN_w
=E[\text{number of row-0 white sites belonging to essential components}].}\tag{3.1}
\]

Horizontal translation symmetry makes the right side

\[
w\rho_w(q).                                                    \tag{3.2}
\]

Hence exactly at every finite width,

\[
\boxed{
\frac{\nu_w}{w}E_WN_w=\rho_w(q).}                            \tag{3.3}
\]

Using (2.3),

\[
\boxed{
\frac{\nu_w}{w}E_WN_w\to\theta_8(q).}                       \tag{3.4}
\]

Combine with (1.5):

\[
\boxed{
\frac{E_WN_w}{wE_WL_w}\to\theta_8(q).}                      \tag{3.5}
\]

This is a ratio-of-component-Palm-means theorem.  It does not yet assert `N_w/(wL_w)->theta` samplewise.

## 4. Boundary mean from the exact dual surface-excess identity

The exact complementary component score on this branch gives

\[
\boxed{
qE_WB_w-pE_WN_w
=pq\,\partial_p\log\nu_w^4(p).}                              \tag{4.1}
\]

The black component-Palm activity estimate gives at every fixed subcritical `p`

\[
|\partial_p\log\nu_w^4(p)|=O_p(w).                            \tag{4.2}
\]

Meanwhile `nu_w` is exponentially small in `w`.  Multiply (4.1) by `nu_w/w`:

\[
q\frac{\nu_wE_WB_w}{w}
-p\frac{\nu_wE_WN_w}{w}
=O_p(\nu_w).                                                   \tag{4.3}
\]

Using (3.4),

\[
\boxed{
\frac{\nu_w}{w}E_WB_w
\to\frac{p}{q}\theta_8(q).}                                 \tag{4.4}
\]

Equivalently,

\[
\boxed{
\frac{E_WB_w}{E_WN_w}\to\frac{p}{1-p}.}                     \tag{4.5}
\]

This recovers the infinite-cluster one-site-flip density

\[
\beta_8(q)=\frac{p}{q}\theta_8(q)                            \tag{4.6}
\]

at the level of giant-cylinder-component Palm means without separately proving a boundary mixing theorem.

## 5. Mean reward vector

Equations (1.5), (3.4), and (4.4) may be collected as

\[
\boxed{
\left(
\nu_wE_WL_w,
\frac{\nu_w}{w}E_WN_w,
\frac{\nu_w}{w}E_WB_w
\right)
\longrightarrow
\left(
1,
\theta_8(q),
\frac{p}{q}\theta_8(q)
\right).}                                                     \tag{5.1}
\]

Thus the three leading Palm means contain no unknown morphology amplitudes.

Since `E L` is order `1/nu`, both additive rewards are order `w/nu`, as expected for an exponentially long supercritical slab.

## 6. Separation from the stronger samplewise bulk law

`supercritical-white-slab-bulk-20260914.md` conjectures

\[
\frac{N_w}{wL_w}\to\theta_8(q),
\qquad
\frac{B_w}{wL_w}\to\frac{p}{q}\theta_8(q)                   \tag{6.1}
\]

in probability under component Palm.

The present theorem is strictly weaker but already nontrivial:

\[
\frac{EN_w}{wEL_w}\to\theta_8(q),
\qquad
\frac{EB_w}{EN_w}\to p/q.                                   \tag{6.2}
\]

Proving (6.1) still requires controlling internal fluctuations/conditioning inside a random exponentially long slab.  It should not be inferred from the mean identities alone.

## 7. Archive-facing checks

A component-Palm simulation of the huge white complementary component can therefore be checked against three levels:

1. already proved: `nu E L ->1`;
2. proved here: `nu E N / w -> theta_8(1-p)`;
3. proved here without knowing `theta`: `E B / E N -> p/(1-p)`.

Only after these mean checks pass should one test the stronger samplewise LLN or Gaussian bulk fluctuations.

## 8. Claim boundary

The Campbell identities are exact.  The local essential-membership limit uses an author-level supercritical white frame built from the standard square-site matching crossing dichotomy plus subcritical black exponential crossing decay.  The full samplewise slab LLN/CLT remains conjectural.