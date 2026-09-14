# Exact 4/8 reflection law for the charge fugacity

2026-09-14.  Short exact consequence of digital Alexander duality, graph inclusion and the charge-neutral coordinates.

## 1. Separate graph coordinates

For graph `G4` or `G8` at occupation density `p`, define

\[
\chi_G(p)=P_G(r=0)+P_G(r=2),                                  \tag{1.1}
\]

\[
\theta_G(p)=\log\frac{P_G(r=2)}{P_G(r=0)}.                   \tag{1.2}
\]

In the rank-one sector let `H_{G,p}(z)` be the conditional essential-component count PGF and let `L_G` denote the projective slope mark.

## 2. Digital Alexander is an exact involution on the full tuple

Configurationwise,

\[
r_4(\omega)+r_8(\omega^c)=2.                                \tag{2.1}
\]

At black density `p`, the complement has matching density `1-p`.  Therefore

\[
P_{8,p}(r=0)=P_{4,1-p}(r=2),                                 \tag{2.2}
\]

\[
P_{8,p}(r=2)=P_{4,1-p}(r=0),                                 \tag{2.3}
\]

and rank one is preserved.

Hence exactly

\[
\boxed{\chi_8(p)=\chi_4(1-p),}                               \tag{2.4}
\]

\[
\boxed{\theta_8(p)=-\theta_4(1-p).}                          \tag{2.5}
\]

The stronger rank-one component theorem gives preservation of count and slope under complement, so

\[
\boxed{H_{8,p}(z)=H_{4,1-p}(z)}                              \tag{2.6}
\]

and the full slope-marked neutral law is transported unchanged.

Thus the exact charge-neutral involution is

\[
\boxed{(\chi,\theta,H,L)_8(p)
=(\chi,-\theta,H,L)_4(1-p).}                                 \tag{2.7}
\]

## 3. Same-occupied-set graph inclusion orders the fugacity

For a fixed occupied set, the matching graph contains the NN graph.  Therefore its ambient rank is at least the NN rank.

At the probability level,

\[
P_{8,p}(r=2)\ge P_{4,p}(r=2),                                \tag{3.1}
\]

\[
P_{8,p}(r=0)\le P_{4,p}(r=0).                                \tag{3.2}
\]

Consequently

\[
\boxed{\theta_8(p)\ge\theta_4(p).}                           \tag{3.3}
\]

For ordinary square tori where the extra matching diagonals produce a strict enhancement with positive probability, the inequality is strict in the interior.

## 4. Exact reflection dominance for the NN charge field

Combine (2.5) and (3.3):

\[
-\theta_4(1-p)\ge\theta_4(p).                                \tag{4.1}
\]

Therefore

\[
\boxed{\theta_4(p)+\theta_4(1-p)\le0.}                       \tag{4.2}
\]

This is the charge-field version of the earlier probability inequalities

\[
P_2(1-p)\le P_0(p),                                          \tag{4.3}
\]

\[
M(p)+M(1-p)\le0.                                              \tag{4.4}
\]

It is stronger conceptually because `theta` is the coordinate whose zero is **exactly** the matching root.

At `p=1/2`,

\[
\boxed{\theta_4(1/2)\le0.}                                  \tag{4.5}
\]

Since `theta_4` is strictly increasing, the finite NN balance root obeys

\[
\boxed{p_*\ge1/2,}                                           \tag{4.6}
\]

with strict inequality under strict matching enhancement.

## 5. Constraint on any scaling crossover

Suppose a family of finite tori admits a charge-field scaling limit after a thermal reparameterization `x=x_w(p)`.  If complement reflection acts asymptotically as `x -> x^c` under `p->1-p`, then the two graph scaling functions must satisfy

\[
\boxed{\Theta_8(x)= -\Theta_4(x^c).}                          \tag{5.1}
\]

Graph inclusion additionally requires the appropriately aligned inequality

\[
\Theta_8(x)\ge\Theta_4(x).                                   \tag{5.2}
\]

A proposed common-window ansatz that violates these two relations is incompatible with finite topology regardless of how well it fits marginal count data.

For a self-matching control model, the involution would reduce to an oddness condition around the self-dual parameter.  The square-site NN/matching pair is not self-matching, so no such extra identification is assumed here.

## 6. Boundary

All identities through Section 4 are exact finite-lattice statements.  Section 5 is only the induced constraint on any future scaling limit.  No critical scaling function is asserted.