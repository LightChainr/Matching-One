# Matching defect polynomial: vertex-subset generating function, not Tutte

Status: C5 identification on axis L=2/3 and C4 N=10. Issue #144 first deliverable. No deletion-contraction algorithm.

## Identification

The finite matching polynomial is not primarily a one-variable root object. Its Bernstein coefficients are

```text
a_k = sum_{|ω|=k} q(ω),
q(ω) = I_wrap(black NN) - I_wrap(white NN+NNN) ∈ {-1,0,+1}.
```

On the implemented square tori, wrapping-difference channels are configuration-identical and the P34 identity

```text
C_black - C_white = q + V - E + F0
```

holds configuration-wise, so `M(p)=E[q]`. The matching polynomial is the Bernoulli generating function of this homology event.

## Tutte obstruction

Do not force the problem into an edge-subset Tutte, Bollobás–Riordan, or Krushkal polynomial.

- Site percolation sums over **vertex subsets**, not edge subsets.
- The black graph is NN and the white graph is NN+NNN, except on a self-matching triangulation.
- Even on C4 N=10, where primal and matching edges coincide, the sum is still over wrapping events of vertex subsets.

A useful negative result is exactly this obstruction, together with the vertex-subset generating function above.

## Exact checks

Axis L=2/3 reproduce the committed matching polynomials from `sum q`. Complement is **not** `q → -q` on the axis (the two graphs differ).

C4 N=10 is self-matching: complement is an involution `q(ω^c)=-q(ω)`, the Bernstein vector is antisymmetric, and the power polynomial is the Beta(3,3) control `12p^5-30p^4+20p^3-1`.

## Boundary

No cheaper deletion-contraction or transfer-matrix algorithm is claimed. Galois complexity of the same polynomials is a separate exact track (#104). Continuum FK/Potts Q-derivatives are a separate theory track (#114).
