# Matching involution and Virasoro-descendant parity

Status: mechanism derivation under an explicit RG assumption.  This note sharpens the parity assignment used in the two-spin-4 model.

## 1. Microscopic involution

For a planar site lattice `G` and its matching lattice `Ghat`, the complementary construction exchanges

\[
(G,p)\longleftrightarrow(\widehat G,1-p).
\]

At the critical pair

\[
p_c(G)+p_c(\widehat G)=1,
\]

both microscopic models flow to the same percolation fixed point.  In a neighborhood of that fixed point the matching/complement operation induces an involution `T` on scaling fields.

The important point is geometric: the operation changes the microscopic connectivity/occupation convention but does not rotate, translate or rescale the continuum coordinates.  The natural RG assumption is therefore

\[
[T,L_n]=[T,\bar L_n]=0,
\]

or, less strongly, that `T` commutes with the local conformal/Euclidean generators acting on a scaling family.

This is an assumption about the scaling-limit action of the matching map, not a theorem of Virasoro algebra alone.  It is testable through the descendant parity spectrum in issue #48.

## 2. Thermal primary is odd

Choose the normalized thermal scaling coordinate `t` so that the two matching models use opposite coordinates:

\[
T:t\mapsto -t.
\]

The continuum perturbation is

\[
S=S_*+t\int\epsilon(x)d^2x+\cdots.
\]

Invariance of the common fixed-point description then requires

\[
\boxed{T\epsilon=-\epsilon}
\]

for the bottom thermal scaling field.  Thus the thermal primary belongs to the matching-odd sector.

## 3. Every ordinary thermal descendant inherits odd parity

If `T` commutes with the Virasoro generators, then for every descendant

\[
|\psi\rangle=L_{-n_1}\cdots L_{-n_k}\bar L_{-m_1}\cdots\bar L_{-m_l}|\epsilon\rangle,
\]

we have

\[
T|\psi\rangle
=L_{-n_1}\cdots\bar L_{-m_l}T|\epsilon\rangle
=-|\psi\rangle.
\]

Therefore the level-4 spin-4 quasiprimary candidate

\[
Q_4^\epsilon\sim
40L_{-2}^2-60L_{-3}L_{-1}-9L_{-4}
\]

and its antiholomorphic conjugate are automatically **matching-odd** if they are physical descendants of the thermal family.

This supplies the missing parity rationale for the `x=21/4` mechanism:

\[
\boxed{x=21/4,\ s=\pm4,\ \eta_T=-1.}
\]

No separate numerical sign assignment is required at the mechanism level; the numerical data test whether the matching map indeed realizes this commuting involution.

## 4. Identity family is even

The identity is fixed:

\[
T\mathbf1=\mathbf1.
\]

Hence ordinary identity-family Virasoro descendants inherit even matching parity.  In particular, any square-lattice spin-4 identity descendant (the analogue of `T^2+\bar T^2` in ordinary CFT language) belongs to

\[
\boxed{\eta_I=+1.}
\]

This gives the structural two-family projector:

```text
identity-family spin 4  -> matching even -> appears at the center in S
thermal-family spin 4   -> matching odd  -> appears at the center in D=M/2
```

The observed finite-size powers may then differ because the first allowed spin-4 descendants have different dimensions.

## 5. Center and derivative selection follows immediately

For a matching-parity eigenfield with `eta`, let its universal thermal dependence be `F(z)`, where matching sends `z->-z`.  Its contributions are

\[
S\propto F(z)+\eta F(-z),
\qquad
D\propto F(z)-\eta F(-z).
\]

Thus

\[
S^{(n)}(0)\neq0\iff(-1)^n=\eta,
\]

\[
D^{(n)}(0)\neq0\iff(-1)^n=-\eta.
\]

For the two candidate families this yields the frozen spectrum already recorded in #48:

```text
identity/even x=4:
  P4[S]   ~ N^-1
  P4[D']  ~ N^-5/8

thermal/odd x=21/4:
  P4[D]   ~ N^-13/8
  P4[S']  ~ N^-5/4
```

The first two center channels are already supported at finite size; the crossed first derivatives are a direct test of the commuting-involution assumption.

## 6. Why this is stronger than an Ising analogy

The Ising irrelevant-operator classification is useful evidence that square lattices can couple to spin-4 descendants in identity and energy families.  But the matching parity assignment here does not need to be copied from Ising duality.

It follows instead from two percolation-specific ingredients:

1. matching/complement reverses the normalized thermal coordinate;
2. the induced scaling-limit map is assumed to commute with conformal descendants.

If #48 violates the derivative parity pattern, the assumption is falsified even if the center exponents remain numerically accurate.

## 7. LCFT/Jordan block: an involution cannot shear a single rank-2 block

At `Q=1` the thermal Kac field belongs to a logarithmic structure.  Logarithmic finite-size amplitudes therefore remain possible, but the involution itself is more constrained than a generic triangular matrix.

Consider one rank-2 Jordan block of `L_0`,

\[
L_0=hI+N,\qquad N^2=0,\quad N\neq0.
\]

If `T` commutes with `L_0`, then on this indecomposable block its commutant has the form

\[
T=aI+bN
\]

(up to basis conventions).  Imposing the exact involution condition

\[
T^2=I
\]

gives

\[
a^2=1,\qquad 2ab=0.
\]

Over characteristic zero, `a=+/-1` and therefore

\[
\boxed{b=0.}
\]

So an involution commuting with `L_0` cannot act as

`partner -> +/- partner + c * bottom`

with nonzero `c` inside a single isolated rank-2 Jordan block.  The two Jordan partners share the same matching parity on that block.

This is important: logarithmic scaling does **not** by itself destroy the matching-parity selection rule.

More complicated mixing can still occur if there are several isomorphic/degenerate indecomposable blocks on which `T` acts nontrivially in the multiplicity space.  That possibility has to respect `T^2=1` and the conformal commutant and should be treated explicitly rather than assumed.

## 8. Why logarithmic finite-size terms are still allowed

Even when `T` acts as `+/- I` on one Jordan block, the non-diagonal action of `L_0` itself produces logarithms under scale evolution.  Therefore a matching-odd thermal logarithmic multiplet can still generate

\[
L^{-13/4}(A+B\log L+\cdots)
\]

while every term remains in the same matching-odd sector.

Consequences:

- parity alternation in #48 is more robust than the no-log amplitude law;
- a nonzero Gaussian-doubling residual can still diagnose a logarithmic/Jordan amplitude;
- such a residual should not be interpreted as matching-parity violation unless the crossed `S/D` derivative rules also fail.

A full LCFT treatment must identify both the indecomposable module and any multiplicity-space action of the matching involution.

## 9. Claim boundary

The following is now a coherent **mechanism deduction**, not yet a theorem:

> If the microscopic matching/complement map induces a coordinate-blind involution of the common percolation fixed point that commutes with Virasoro generators, then identity-family descendants are matching-even and thermal-family descendants are matching-odd.  The proposed x=21/4 spin-4 thermal quasiprimary therefore has precisely the parity required to control the central matching residual.  Within a single LCFT Jordan block, the same assumptions force the whole block to carry one matching parity rather than a triangularly sheared action.

What would promote this further:

- an explicit FK/Temperley-Lieb realization of `T` in the scaling state space;
- the derivative-parity predictions of #48 passing prospectively;
- self-matching/self-dual controls (#42/#44) showing the expected parity projection;
- a consistent identification of the Q=1 indecomposable module and any multiplicity-space action.
