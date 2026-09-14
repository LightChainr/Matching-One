# Canonical birth-gap amplitude: a thermal-metric-free scalar from the whole rank curve

Date: 2026-09-14. Exact identity plus finite high-precision diagnostics for draft PR #773.

Let `T1<=T2` be the two homology-rank birth times in the monotone coupling. Let

    b(p)=1/2 log[P0(p)/P2(p)]

be the exact canonical rank-odds coordinate. It is strictly decreasing from `+infinity` to `-infinity` along the homogeneous Bernoulli family.

Define the **canonical birth gap**

    G_b = b(T1)-b(T2) >= 0.                               (1)

Unlike the microscopic gap `T2-T1`, this quantity is invariant under arbitrary monotone reparameterizations of p because b itself is defined by the rank law.

## 1. Exact one-time formula

For each realization,

    b(T1)-b(T2)
      = int_(T1)^(T2) [-b'(p)] dp.                        (2)

Using `P1(p)=P(T1<=p<T2)` and Tonelli,

    E G_b
      = int_0^1 [-b'(p)] P1(p) dp                         (3)

and, changing variables from p to b,

    E G_b
      = int_(-infinity)^(infinity) P1(b) db.              (4)

In the softmax coordinates `(b,d)`,

    P1(b)=e^{d(b)}/[e^b+e^{d(b)}+e^{-b}],                 (5)

so the canonical mean gap is a scalar functional of the **entire self-normalized rank curve** `d=D_L(b)`.

It therefore needs no pc, no thermal metric and no derivative normalization.

## 2. Source-scan interpretation

The topological-source balance scan satisfies `b=s`. Hence the same quantity is

    E G_b = int_R P1,balance(s) ds,                        (6)

where `P1,balance(s)` is the rank-one probability after applying source s and retuning p so the sourced rank-0 and rank-2 masses are equal.

Thus an all-width closure/source formalism can compute the canonical gap directly by integrating a one-parameter family of sourced balance problems.

## 3. Finite exact-polynomial diagnostics

Using the committed exact rank-by-occupation polynomials and high-precision quadrature of the exact rational integrand `[-b'(p)]P1(p)`, the current controls are

| lattice | L=3 | L=4 |
|---|---:|---:|
| square NN / matching rank law | 1.1119239786199567 | 1.1730803897060144 |
| triangular self-matching | 1.1494745051701594 | 1.1803623670395626 |

These values are numerical evaluations of exact finite polynomials, not interval-certified integrals. The L=4 cross-lattice difference is about 0.6%. Two widths do **not** prove universality, but the agreement is strong enough to justify making this a frozen output in #775/#776.

For comparison, the microscopic p-coordinate mean gaps are exactly

| lattice | L=3 | L=4 |
|---|---:|---:|
| square | `3/20 = 0.15` | `14122/109395 = 0.1290918232...` |
| triangular | `3/20 = 0.15` | `991/7735 = 0.1281189399...` |

The microscopic gaps require a thermal metric to compare across lattices; the canonical gap (4) does not.

## 4. Continuum target

If the self-normalized rank curve converges at fixed torus modulus,

    D_L(b) -> D_*(b;tau),                                 (7)

and dominated convergence can be justified in the tails, then

    G_*(tau)
      = int_R
          exp[D_*(b;tau)]
          / [2 cosh b + exp D_*(b;tau)]
        db                                                   (8)

is a universal topological scalar for that modulus.

This quantity is not determined by the single critical wrapping probability `a_*(tau)`; it probes the full near-critical rank equation of state. Hence it is complementary to the Pinson critical-point control.

Self-matching triangular site is the cleanest rigorous/numerical positive-control route because its finite curve is exactly even.

## 5. Higher canonical gap moments

The one-time rank curve determines only the first canonical gap moment. For k>=2,

    E[G_b^k]
      = int_{R^k}
          P( all queried b_i lie between b(T2) and b(T1) )
        db_1...db_k,                                      (9)

which requires multi-source/multi-time rank-one persistence, just as higher moments of `T2-T1` require the two-time monotone coupling surface.

No independence of the two births is assumed.

## 6. Recommended output

Exact rank-sector transfer tasks #775 and triangular control #776 should report:

1. the microscopic exact mean gap `int P1(p)dp`;
2. the canonical mean gap (3)/(4) at high precision;
3. the source-balanced integrand `P1,balance(s)` on the same frozen s-grid used for the self-normalized odd residual.

The canonical gap should be treated as a **parameterization-free universality diagnostic**, not fitted as a new critical exponent.
