# A self-normalized rank equation of state: eliminate the thermal coordinate entirely

Date: 2026-09-14. Exact finite construction plus a universality target for draft PR #773.

Start from the canonical rank-simplex coordinates

    b(p) = (1/2) log[P_0(p)/P_2(p)],
    c(p) = P_1(p)/[2 sqrt(P_0(p)P_2(p))].                 (1)

No continuum assumption is needed to define them.

## 1. `b` is a global finite thermal/topological coordinate

Use logit `z=log[p/(1-p)]`. Sector thermodynamics gives

    db/dz = -(1/2)[E(K|r=2)-E(K|r=0)].                    (2)

The rank-2 event is increasing and the rank-0 event decreasing under site addition. Since K is increasing, Harris association gives

    E(K|r=2) >= E K >= E(K|r=0),                          (3)

with strict inequalities for a nondegenerate honest torus and `0<p<1`. Therefore

    db/dz < 0.                                            (4)

As `p->0`, `P_0->1` and `P_2->0`, so `b->+infinity`; as `p->1`, `b->-infinity`. Hence `b` is a one-to-one global coordinate along the physical homogeneous Bernoulli family.

This is stronger than merely saying the balance root is unique. Every occupation probability can be reparameterized by the source required to rebalance its extreme topology sectors.

## 2. Eliminate `p`, `lambda`, and the thermal metric

Because b is invertible, define the exact finite curve

    c = C_L(b).                                           (5)

This curve is invariant under **any** monotone reparameterisation of the microscopic occupation coordinate, including nonlinear thermal normal-coordinate choices. It therefore separates two questions cleanly:

- locating `p_c` and the finite balance root requires the p-to-thermal map;
- the shape of the rank law once parametrized by its own odds coordinate b does not.

For the primal/matching pair,

    b_G(p) = - b_Ghat(1-p),
    c_G(p) =   c_Ghat(1-p),                               (6)

so their self-normalized curves are mirrors. If both microscopic models converge to one continuum near-critical rank law, the universal curve must satisfy

    C_*(b)=C_*(-b).                                       (7)

Thus odd parts of the finite `C_L(b)` are a direct measure of matching-asymmetric corrections.

## 3. Parameterization-free Taylor invariants at balance

At the balance point b=0, define conditional occupation cumulant combinations

    g_n = kappa_n(K|2)-kappa_n(K|0),
    h_n = kappa_n(K|1)
          -(1/2)[kappa_n(K|0)+kappa_n(K|2)].               (8)

Then

    b_z = -g_1/2,
    b_zz = -g_2/2,
    (log c)_z = h_1,
    (log c)_zz = h_2.                                    (9)

Chain rule gives the exact first derivative

    d(log c)/db = -2 h_1/g_1,                             (10)

and exact second derivative

    d^2(log c)/db^2
      = 4 (g_1 h_2 - g_2 h_1)/g_1^3.                     (11)

Unlike `g_2/g_1`, (11) is invariant not only under an overall thermal metric but under arbitrary smooth reparameterization of p: it is an intrinsic curvature of the rank-simplex curve.

Exact controls from the committed L=3,4 rank histograms are

| L | `d log c/db` | `d^2 log c/db^2` |
|---|---:|---:|
| 3 | 0.06240221458 | -0.12565256793 |
| 4 | 0.03379345107 | -0.11225277224 |

Two sizes establish no limit. The decreasing first derivative is merely compatible with the expected even continuum curve. The second derivative is a candidate universal shape number for the square-torus rank law.

## 4. Relation to the near-critical Taylor jet

If a conventional thermal coordinate lambda is used, write

    b(lambda)=b_1 lambda + b_3 lambda^3/6 + ...,
    ell(lambda)=log c(lambda)=ell_0+ell_2 lambda^2/2+....   (12)

Then (11) tends simply to

    ell_2/b_1^2.                                          (13)

So the intrinsic curve curvature is the metric-free even Taylor invariant previously written using derivatives with respect to z. It can be compared across lattices without first solving the mass-normalization problem in #770.

By contrast, the odd invariant `b_3/b_1^3` still depends on the choice of nonlinear thermal coordinate, although it is invariant under an overall linear metric factor. This distinction is important: `C(b)` is the cleaner object whenever only topology-sector probabilities are being compared.

## 5. Matching-odd irrelevant-field cross-test

Suppose, as a conditional one-field model, the finite theory contains one leading matching-odd irrelevant amplitude

    u_L ~ L^(-omega),

and the continuum thermal exponent is `y_t=3/4`. Symmetry permits the expansion

    b(lambda,u)=B(lambda)+u H(lambda)+...,
    log c(lambda,u)=C(lambda)+u J(lambda)+...,             (14)

where `B,J` are odd and `C,H` even. Then

    b(0,u)=u H(0)+...,

so the root displacement in microscopic p has exponent

    omega_root = omega + y_t.                              (15)

At the shifted balance root the even-coordinate slope scales as

    (log c)_z = O(L^(-omega+y_t)),                         (16)

hence

    omega_root - omega_c-slope = 2 y_t = 3/2.             (17)

This gives an independent multi-observable test of any proposed leading matching-odd correction.

Examples:

- 8-arm proposal: root `L^-4` -> `omega=13/4` -> `(log c)_z ~ L^-5/2`;
- 6-arm proposal: root `L^-5/3` -> `omega=11/12` -> `(log c)_z ~ L^-1/6`.

Multiple irrelevant fields, vanishing amplitudes or logarithmic partners can spoil a one-power fit; (17) is therefore a conditional consistency relation, not a theorem about the actual lattice.

The exact finite values

    (log c)_z = -0.09188772912...  (L=3),
                 -0.06302555101...  (L=4)

are recorded only as future controls. Issue #769 has been asked to output the L=5 value at negligible additional cost.

## 6. Use in continuum/map-resolved work

For each torus modulus/twist tau, `C_*(b;tau)` is a fully self-normalized topological equation of state. It is a stronger positive control than one critical probability because it tests a whole local curve while keeping explicit rank/connectivity semantics.

This construction does not identify a CFT field or solve original-U. It provides a parameterization-free target that a proposed map-resolved continuum basis should reproduce before more delicate source amplitudes are interpreted.
