# Topological Clapeyron response: root motion as a difference of quasi-stationary source scores

Date: 2026-09-14

Status: exact finite-state Perron/Feynman--Hellmann calculus once a physical source is represented in the two safe transfer operators.  Large-width/CFT interpretations are separate.

## 1. The matching root is a coexistence line between two topological void phases

At fixed circumference `w`, let

```text
I4(p,g) = -log lambda4(p,g),
I8(q,g) = -log lambda8(q,g),
q=1-p,
```

be the Perron free energies of the NN-safe and complementary-matching-safe transfer operators in the presence of a microscopic source `g`.

The source must be physically specified on each graph/colour representation; the symbols `g` on the two sides mean the declared matched perturbation, not an arbitrary matrix entry.

Define the charge free-energy difference

```text
Theta_w(p,g)=I4(p,g)-I8(1-p,g).
```

The fixed-width charge-coexistence root `p_w(g)` satisfies

```text
Theta_w(p_w(g),g)=0.
```

Thus the matching root is literally a coexistence line between two quasi-stationary topological void phases.

## 2. Perron/Doob source scores

Let a safe transfer entry carry microscopic weight `W_G(i,M,j;p,g)`.  Normalize positive Perron vectors by `l^T r=1` and define the one-step Doob/Perron law

```text
Q_G(i,M,j)
 = l_i W_G(i,M,j) r_j / lambda_G.
```

For a source derivative define the physical row score

```text
H_G = partial_g log W_G |_(g=0).
```

Feynman--Hellmann gives exactly

```text
partial_g log lambda_G = E_QG[H_G],
partial_g I_G = -E_QG[H_G].
```

Therefore

```text
Theta_g
 = -E_Q4[H_4] + E_Q8[H_8].
```

No rare torus probability needs to be estimated.

## 3. Exact topological Clapeyron equation

Implicit differentiation of the coexistence equation gives

```text
boxed:
dp_w/dg
 = -Theta_g/Theta_p
 = [E_Q4 H_4 - E_Q8 H_8]/Theta_p.
```

The thermal denominator is already known exactly from the row-occupation Perron identity:

```text
Theta_p
 = [w-Kbar_4^0(p)-Kbar_8^0(1-p)]/[p(1-p)] > 0.
```

Hence the root moves if and only if the physical source has a nonzero expectation difference between the two safe quasi-stationary phases.

This is the finite-width version of the question posed in the research compass:

> which microscopic correction is genuinely common to the two topological sectors, and which one has a nonzero difference matrix element?

The answer can be measured before naming a continuum field.

## 4. Common versus differential source directions

For every source define its two phase scores

```text
mu_4(g)=E_Q4 H_4,
mu_8(g)=E_Q8 H_8.
```

Decompose them into

```text
mu_common = (mu_4+mu_8)/2,
mu_diff   = (mu_4-mu_8)/2.
```

Only `mu_diff` moves the coexistence root at first order:

```text
dp_w/dg = 2 mu_diff / Theta_p.
```

Thus a large correction in each individual sector can be irrelevant to Matching One if it lies almost entirely in the common direction.

This gives an operational meaning to the observed hierarchy

```text
common x≈4 correction >> sector-odd mismatch.
```

The common `x≈4` field can be numerically large in `I4` and `I8` yet cancel from `Theta`.  Its mixed dressing of a differential spin-four source can re-enter at the next order, as proposed in `sector-even-dressing-of-spin4-tower-20260914.md`.

## 5. The oblique spin-four result is already a Clapeyron numerator measurement

Changing the orientation of the cylinder relative to the square lattice changes the matrix element of an anisotropic lattice correction while preserving the thermal parameter.

At criticality, the existing oblique safe-transfer data show

```text
Theta_w(pc,theta)
 ~ B4 cos(4theta) ell^(-17/4)
```

with very small angular-even leakage.

In the Clapeyron language this says that the effective sector-difference numerator is dominated by one spin-four angular direction.  The root law

```text
p_root-pc ~ -A4 cos(4theta) ell^-4
```

then follows after division by the leading scalar thermal denominator.

The exact-equal-circumference `(4,3),n=2` versus axis `w=10` experiment therefore tests numerator and denominator separately:

```text
E_43/E_axis ~ cos4theta_43,
D_43/D_axis ~ 1,
root-shift ratio ~ cos4theta_43.
```

A failure can be localized to either source numerator or thermal denominator rather than described vaguely as a bad exponent fit.

## 6. Relation to the `(T,N)` response of #802

The root tangent in the canonical rank coordinates is

```text
T_g = dp*/dg = -b_g/b_p.
```

The Clapeyron formula is the fixed-width safe-transfer realization of the same tangent when `b` is represented by the topological free-energy difference.

For a second shape observable `e`, the fixed-b normal response is

```text
N_g=e_g-(e_p/b_p)b_g.
```

Thus an actual source should be summarized by the pair

```text
(T_g,N_g).
```

The new point is that `T_g` itself can be decomposed into a **difference of phase scores** before any continuum extrapolation.  This makes the first coordinate much more interpretable.

A source can have:

```text
mu_diff !=0 but N_g=0 : pure coexistence-line/thermal-tangent motion;
mu_diff =0 but N_g!=0 : shape change without first-order root motion;
both nonzero           : genuine mixed response.
```

These cases should not be conflated.

## 7. A direct route for real spatial or local sources

For a site-local/logit perturbation the row score `H_G` is explicit.  For example, a logit field on newly exposed sites has score equal to the corresponding occupied count minus its normalization contribution; source-compatible safe transfer therefore gives `mu_4,mu_8` directly by one Perron left/right contraction.

The finite spatial-source Hessian work already shows that mean-zero sources often have zero first-order response by translation symmetry.  In that case the first nonzero coexistence response is second order and should be computed from the differentiated Perron resolvent / Green--Kubo susceptibility, not by reinterpreting a symmetry zero as a field selection rule.

Thus the source order (linear versus quadratic) must be fixed before comparing spin channels.

## 8. Second-order coexistence curvature

If symmetry forces `Theta_g=0`, differentiating twice gives at `g=0`

```text
p_w''(0) = -Theta_gg/Theta_p
```

when `p_w'(0)=0`.

For a transfer source with no explicit second derivative in the log weight, `Theta_gg` is the difference of integrated score covariances (Green--Kubo susceptibilities) between the two Doob phases, plus any declared contact term from the physical coordinate.

This is the semi-infinite-cylinder analogue of the finite-torus conditional-covariance root Hessian already derived in the spatial-source work.

It suggests a useful division of labour:

```text
finite torus  : rank-conditioned covariance difference;
long cylinder : safe-phase Green--Kubo susceptibility difference.
```

A continuum source identification should make these compatible in their common scaling regime.

## 9. Topological rank source and Krushkal channel

For bond/FK surface states, the rank charge `h` has the exact Krushkal realization

```text
exp[h(r-1)] = A^(s/2) B^(s_perp/2),
A=e^h, B=e^-h.
```

Therefore the differential phase score can in principle be formulated as a mixed response

```text
partial_g partial_h F_topological.
```

This gives a source-defined version of the phrase “sector-odd matrix element”.  A candidate continuum field matters for the matching root only if it contributes to this mixed topological response.

## 10. A high-information source matrix rather than another exponent table

Suppose we choose a small declared basis of real microscopic perturbations `g_a`:

```text
uniform thermal,
spin-four anisotropy,
one scalar local perturbation,
one map/connectivity-resolved perturbation when available.
```

For each width compute only

```text
D_a(w) = E_Q4 H_a - E_Q8 H_a,
C(w)   = Theta_p.
```

The vector `D_a` is the finite microscopic **differential-source fingerprint** of the two sectors.  Its angular/transformation properties can be compared across widths and geometries.

This is more informative than fitting one root sequence to several possible powers because it asks directly which physical source directions survive the common-sector quotient.

No large source dictionary is justified: two mechanisms should be selected first, and one perturbation should be chosen that gives different predicted fingerprints.

## 11. Implication for original-U

This does not solve #275 or replace its frozen observable.  It does clarify the theoretical object required there.

A named candidate must eventually provide a forward map from the actual original-U microscopic source into a differential topological response, including normalization and moving-root counterterm.  The finite-width Clapeyron score is a lattice-side prototype of that map:

```text
physical source
 -> phase-score difference
 -> root tangent.
```

If two continuum candidates give the same allowed differential-source image, more precision on the same source cannot identify them.

## 12. Claim boundary

Exact finite-state statements:

- Perron/Doob source score identity;
- coexistence implicit derivative;
- thermal denominator from row occupation;
- common/differential decomposition.

Conjectural/scaling interpretations:

- identification of particular differential source directions with spin-four/eight-arm continuum operators;
- large-width scaling of the score vector;
- compatibility with original-U continuum forward maps.

The main conceptual result is independent of those identifications: **Matching One responds to the difference of source scores between two topological quasi-stationary phases.**