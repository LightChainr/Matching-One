# Parameterization-free matching-odd residual: test the irrelevant exponent without pc or a thermal metric

Date: 2026-09-14. Additive analysis note for draft PR #773.

The finite rank law has exact canonical coordinates

    b = (1/2) log(P0/P2),
    c = P1/[2 sqrt(P0 P2)].

Since b is strictly monotone along the homogeneous Bernoulli family, every finite lattice defines an intrinsic rank curve

    c = C_L(b),                                             (1)

with no reference to p, pc, a thermal metric, or a fitted scaling variable.

## 1. Exact matching-even / matching-odd decomposition of the curve

Digital Alexander complement gives, for the primal/matching pair,

    C_G,L(b) = C_Ghat,L(-b).                                (2)

Therefore the intrinsic logarithmic curve has the exact pair decomposition

    E_L(b) = 1/2 [log C_L(b)+log C_L(-b)],
    O_L(b) = 1/2 [log C_L(b)-log C_L(-b)].                 (3)

`E_L` is even and `O_L` odd by construction. In a self-matching lattice such as the standard triangular-site control, `O_L(b)=0` **exactly at every finite size**.

For square site, `O_L` is a pure finite matching-asymmetry observable. It is not contaminated by a nonlinear choice of the microscopic thermal coordinate, because p has already been eliminated in favour of b.

At b=0,

    O_L'(0) = d(log c)/db
            = -2 h1/g1,                                    (4)

where

    g1 = E[K|2]-E[K|0],
    h1 = E[K|1]-(E[K|0]+E[K|2])/2.

The L=3,4 square exact controls are

    O_3'(0)=0.06240221458...,
    O_4'(0)=0.03379345107....                               (5)

Two points do not determine an exponent.

## 2. Direct relation to a matching-odd irrelevant field

Assume, conditionally, that after using the intrinsic b coordinate there is one leading matching-odd irrelevant scaling field with amplitude

    u_L ~ L^(-omega).                                      (6)

The most general local form compatible with complement parity is

    log C_L(b)
      = F_even(b)
        + u_L J_odd(b)
        + higher terms,                                    (7)

with F_even(-b)=F_even(b) and J_odd(-b)=-J_odd(b). Then

    O_L(b) = u_L J_odd(b)+o(u_L).                          (8)

Thus **the intrinsic odd residual measures omega directly**. There is no additional thermal exponent y_t.

By contrast, the finite balance-root displacement in p has

    p_L^*-pc ~ L^[-(omega+y_t)]                            (9)

under the same one-field model. Therefore

    root exponent - intrinsic-odd exponent = y_t=3/4.     (10)

The earlier z-coordinate slope `h1=(log c)_z` has exponent `omega-y_t`; the b-coordinate residual is cleaner because the thermal Jacobian is absent.

## 3. 6-arm versus 8-arm predictions

For the proposed 8-arm matching-odd mechanism,

    root exponent = alpha_8-alpha_4 = 4,
    omega_8 = 4-y_t = 13/4.                                (11)

Hence

    O_L(b) ~ L^(-13/4).                                    (12)

For a surviving 6-arm mechanism,

    root exponent = alpha_6-alpha_4 = 5/3,
    omega_6 = 5/3-y_t = 11/12,                             (13)

hence

    O_L(b) ~ L^(-11/12).                                   (14)

These powers are so different that widths 5--8 may be able to reject one scenario without pretending to measure a precise asymptotic exponent. Multiple odd irrelevant fields, logarithmic partners, or a vanishing leading amplitude remain explicit alternatives.

## 4. Recommended exact-transfer diagnostic

Issue #775 is already tasked with exact rank-sector occupation polynomials through bounded widths. From each polynomial table, construct C_L(b) by using b itself as the interpolation coordinate. Report O_L(b) at a **predeclared fixed b grid**, for example

    b in {0.1,0.2,0.3,0.4},                               (15)

provided all widths reach those b values within a numerically safe p interval. The grid must be frozen before seeing width trends.

For each b report

    O_L(b),
    E_L(b),
    O_L(b)/b,
    O_L'(0) from exact sector moments.                     (16)

Do not fit a free multi-exponent model. A useful incompatibility screen is whether ratios across widths are even remotely consistent with L^-13/4 versus L^-11/12.

The triangular self-matching control #776 must return O_L(b)=0 to arithmetic precision for every tested b. This is a stronger regression than checking only the balance point.

## 5. Why this is preferable to another pc fit

The observable (3):

- needs no numerical pc;
- needs no mass metric or nonlinear p->lambda map;
- is exactly typed by matching complement parity;
- uses the full rank triplet rather than only the root of P2-P0;
- separates the correction exponent omega from the thermal exponent y_t.

A failure of the 8-arm power in O_L does not by itself identify the replacement field, but it would directly invalidate the simplest one-field 8-arm explanation of the L^-4 root law.

No claim in this note upgrades the 8-arm hypothesis before #768's map/selection audit.
