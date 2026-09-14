# The self-normalized rank curve is exactly a topological-source balance scan

Date: 2026-09-14. Exact finite identity for draft PR #773.

Let the unsourced homogeneous rank probabilities at occupation parameter p be `P_j(p)`, with canonical coordinates

    b(p)=1/2 log[P0(p)/P2(p)],
    c(p)=P1(p)/[2 sqrt(P0(p)P2(p))].                       (1)

Couple a finite topological source s to `X=r-1`:

    P0 -> P0 e^(-s)/Z_s,
    P1 -> P1/Z_s,
    P2 -> P2 e^(s)/Z_s.                                   (2)

As shown in the canonical-coordinate note,

    b -> b-s,
    c -> c.                                                (3)

## 1. Source-balanced scan theorem

For each real source s, define `p_*(s)` as the homogeneous occupation probability at which the **sourced** rank law is balanced,

    P0(p)e^(-s) = P2(p)e^(s).                              (4)

Since b(p) is strictly monotone, this root is unique and (4) is exactly

    b(p_*(s))=s.                                           (5)

The even coordinate is unchanged by the source. Hence

    c_sourced,balance(s)
       = c(p_*(s))
       = C_L(s),                                           (6)

where `c=C_L(b)` is the parameterization-free intrinsic rank curve of the unsourced Bernoulli family.

> **Exact conclusion.** Scanning a topological source and rebalancing the thermal parameter directly traces the intrinsic curve `C_L(b)`, with the source value itself equal to the canonical coordinate b. No numerical inversion of p, knowledge of pc, or thermal metric is needed.

## 2. Direct sourced form of the matching-odd residual

The parameterization-free odd residual becomes

    O_L(s)=1/2 [log c_bal(s)-log c_bal(-s)].                (7)

Thus the proposed irrelevant-exponent diagnostic can be implemented as two sourced balance calculations at `+s` and `-s`.

For a self-matching lattice, complement symmetry gives

    c_bal(s)=c_bal(-s)                                     (8)

at every finite size, so `O_L(s)=0` exactly. This is an especially strong source-pipeline regression for triangular site.

For square site, under the one-leading-matching-odd-field model,

    O_L(s) ~ L^(-omega) J(s).                              (9)

The 8-arm and 6-arm candidates predict respectively `omega=13/4` and `11/12`.

## 3. Transfer-matrix interpretation

The source s does not alter the local site law or local transfer operator. It only changes topology-sector closure weights:

    Z_s = e^(-s) Z_0D + Z_1D + e^s Z_2D.                 (10)

Therefore an all-width closure formalism can compute the source scan without enlarging the local propagation state merely to encode s. This is a natural target for #636 if that reserve line is activated.

At each source s, solve the single scalar balance condition

    e^s Z_2D(p) - e^(-s) Z_0D(p) = 0,                    (11)

then evaluate

    c_bal(s)= Z_1D(p_*)/[2 sqrt(Z_0D(p_*) Z_2D(p_*))].    (12)

The common Bernoulli normalizer cancels.

This construction distinguishes closure-source response from a microscopic field insertion. Repeated zeros or source derivatives obtained from (10) are therefore not by themselves local-operator Jordan diagnostics.

## 4. Derivatives of the sourced balance curve

Because `p_*(s)` is the inverse of b(p), its logit-coordinate derivatives are the inverse-function formulas already derived from conditional occupation cumulants. In logit z,

    dz_*/ds = -2/g1,
    d2z_*/ds2 = -4g2/g1^3,
    d3z_*/ds3 = 8(g1 g3-3g2^2)/g1^5,                     (13)

where `g_n=kappa_n(K|2)-kappa_n(K|0)` at the unsourced balance point.

The derivative of the **shape** along the sourced balance scan is simpler:

    d log c_bal/ds = d log C_L/db,                         (14)

so at s=0 it is exactly the intrinsic matching-odd slope `-2h1/g1`. Higher source derivatives give the intrinsic curve Taylor coefficients without thermal-coordinate contamination.

## 5. Practical recommendation

For exact transfer results in #775, retain both implementations:

1. reconstruct `C_L(b)` from rank-sector polynomials as a function of p;
2. independently solve the sourced balance equation (11) at a frozen source grid and evaluate (12).

They must agree. This is an internal check that simultaneously tests rank-sector normalization, source weights, root finding and the self-normalized coordinate transformation.

A suggested frozen grid is `s in {0, +/-0.1, +/-0.2, +/-0.3, +/-0.4}` if all widths remain numerically safe. The exact grid should be fixed before inspecting width trends.

No asymptotic claim follows from the identity itself.
