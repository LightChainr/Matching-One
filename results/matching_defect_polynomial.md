# Matching defect polynomial (vertex-subset, not Tutte)

Source: `scripts/matching_defect_polynomial.py`.
Claim level: C5 identification `M=E[q]` on the enumerated quotients; C5 obstruction
against an edge-subset Tutte specialization. Issue #144 first deliverable.

## Identification

```text
q(ω) = I_wrap(black NN) - I_wrap(white NN+NNN) ∈ {-1,0,+1}
a_k  = sum_{|ω|=k} q(ω)
M(p) = sum_k a_k p^k (1-p)^{N-k} = E[q]
```

On square tori this is the same `q` as in the P34 identity
`C_black-C_white = q+V-E+F0`.

## Axis enumerations

### axis L=2, N=4

Bernstein `a_k`: `[-1, -4, -2, 4, 1]`
Power basis: `-2*p^4 + 4*p^2 - 1`
P34 identity failures: 0
Complement involution failures: 2
Primal equals matching: False

### axis L=3, N=9

Bernstein `a_k`: `[-1, -9, -36, -78, -90, -36, 36, 36, 9, 1]`
Power basis: `-4*p^9 + 18*p^8 - 18*p^7 + 6*p^3 - 1`
P34 identity failures: 0
Complement involution failures: 246
Primal equals matching: False

## C4 self-matching N=10

Primal and matching graphs coincide. Occupation complement is an involution
`q(ω^c)=-q(ω)`, so `M(1-p)=-M(p)`.

Bernstein `a_k`: `[-1, -10, -45, -100, -100, 0, 100, 100, 45, 10, 1]`
Power basis: `12*p^5 - 30*p^4 + 20*p^3 - 1`
Complement involution failures: 0

## Tutte obstruction

Site matching is a vertex-subset model. The black and white graphs are NN
versus NN+NNN except on a self-matching triangulation, and even there the
sum is over wrapping events of vertex subsets, not over edge subsets of a
single ribbon graph. No Bollobás–Riordan / Krushkal specialization is used.

A cheaper deletion-contraction algorithm is **not** established.
