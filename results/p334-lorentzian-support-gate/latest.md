# Exact gates on the projective-line ULC proof routes

## 1. The multiaffine Lorentzian route closes at the smallest quotient

At N=4, matrix `[[2, 0], [0, 2]]`, fixed line `[0, 1]`, the family is `[3, 12]`.
The two support exponents `[2, 1, 1, 0, 0]` and `[2, 0, 0, 1, 1]` fail symmetric exchange at coordinate 1; every allowed target leaves at least one exchanged exponent outside the fixed-line support.
Because M-convex support is necessary for a homogeneous Lorentzian polynomial with nonnegative coefficients, the proposed homogenized multiaffine polynomial is not Lorentzian in general. No Hessian test is needed after this support failure.

## 2. Natural rank-polynomial real-rooted strengthenings also fail

The normalized-q generating polynomial first fails at N=6, `[[2, 0], [0, 3]]`, line `[1, 0]`. After its zero factor is removed its coefficients are `['1/5', '3/5', '3/5']` and it has 0 real roots out of squarefree degree 2.
The raw count generating polynomial first fails at N=8, `[[2, 0], [0, 4]]`, line `[1, 0]`. Exact Sturm variations leave 2 real roots out of squarefree degree 4.
Thus the observed ULC cannot be promoted to either of these standard real-rooted statements.

## 3. Even normalized matching is too strong

The first exact normalized-matching failure occurs at N=11, `[[11, 3], [0, 1]]`, line `[1, -3]`, between layers 7 and 8.
The maximum flow is 605/726; the violating cut is `11/66 > 0/11`. In particular mask 471 is a same-line rank-one state with no same-line one-site extension even though the next layer contains 11 states; each missing-site insertion jumps to rank two.

## Mechanism classification

- **Exact closure:** multiaffine Lorentzian support, normalized-q real-rootedness, raw-count real-rootedness, and full normalized matching are all false in this family.
- **Still exact bounded evidence:** the normalized layer sequence remains ULC for every line checked through N=12 in the preceding atlas.
- **Revised conjecture:** fixed-line ULC is a rank-sequence phenomenon weaker than the standard support and layer-expansion certificates above.
- **Next proof target:** derive a direct two-step injection or coefficient inequality for `A_k^2 binom(N,k-1) binom(N,k+1) >= A_(k-1) A_(k+1) binom(N,k)^2`; do not seek a global exchange axiom.
