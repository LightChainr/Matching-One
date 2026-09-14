# The safe/void transfer as the periodic dilute-TL magnetic `alpha=0` sector

2026-09-14.  Literature-grounded sector dictionary plus independent square-site transfer checks.  This note sharpens `safe-frontier-dilute-tl-bridge-20260914.md`: the common connectivity module is not only a central-trinomial zero-string space.  The **particular twist selected by the no-horizontal-homology condition** is naturally

\[
\boxed{d=0,\qquad \omega=e^{i\gamma}=i,\qquad \alpha=\omega+\omega^{-1}=0.} \tag{1}
\]

At criticality this standard module has exactly the percolation magnetic conformal weight.  For triangular-site percolation the periodic dilute Temperley--Lieb representation is an integrable lattice theorem.  For the square-site Bernoulli safe transfer the claim is a topological/module dictionary, not local Yang--Baxter integrability.

## 1. The published periodic dilute-TL module

Morin-Duchesne, Kluemper and Pearce study critical triangular-site percolation through the Yang--Baxter-solvable dilute `A_2^(2)` loop model.  Its periodic standard modules are labelled by

\[
W_{N,d,\omega},\qquad \omega=e^{i\gamma},                    \tag{1.1}
\]

where `d` is the number of defects/through-lines and `omega` is the periodic twist.

For the periodic ground state they obtain

\[
(h,\bar h)
=
\left(
\Delta_{\gamma/\pi,d/2},
\Delta_{\gamma/\pi,-d/2}
\right),                                                     \tag{1.2}
\]

with

\[
\Delta_{r,s}
=\frac{(3r-2s)^2-1}{24}.                                    \tag{1.3}
\]

In the zero-defect module the noncontractible loop fugacity is

\[
\boxed{\alpha=\omega+\omega^{-1}=2\cos\gamma.}               \tag{1.4}
\]

The same paper's zero-defect dimensions are the central trinomial numbers.

Primary source: A. Morin-Duchesne, A. Kluemper, P. A. Pearce, *Critical site percolation on the triangular lattice: From integrability to conformal partition functions*, arXiv:2211.12379v2.

## 2. Why the safe transfer chooses `alpha=0`

The square-site safe transfer retains a frontier state until an occupied component closes a horizontally noncontractible cycle.  The moment such a cycle appears the transition is rejected.

In a periodic loop representation this is exactly the specialization in which a noncontractible closed loop has zero statistical weight:

\[
\boxed{\alpha=0.}                                             \tag{2.1}
\]

Equation (1.4) then gives

\[
\omega=\pm i,
\qquad
\gamma=\pm\frac\pi2\pmod{2\pi}.                             \tag{2.2}
\]

The two signs are conjugate descriptions of the same real zero-loop-weight condition.  Choose `omega=i`, `gamma=pi/2`.

The frontier itself carries no propagating noncontractible defect/string, so

\[
\boxed{d=0.}                                                   \tag{2.3}
\]

Thus the topological safe condition picks the published module

\[
\boxed{W_{N,0,i}.}                                             \tag{2.4}
\]

This is much sharper than saying only that the state count resembles a dilute zero-string module.

## 3. The magnetic exponent follows without fitting

Set

\[
r=\gamma/\pi=1/2,
\qquad s=d/2=0.                                               \tag{3.1}
\]

Then

\[
\Delta_{1/2,0}
=\frac{(3/2)^2-1}{24}
=\frac5{96}.                                                  \tag{3.2}
\]

Therefore

\[
\boxed{h=\bar h=5/96,
\qquad x=h+\bar h=5/48.}                                     \tag{3.3}
\]

This is precisely the percolation magnetic scaling dimension observed independently in the safe Perron excitation:

\[
w I^0_w(p_c)\longrightarrow2\pi\frac5{48}.                  \tag{3.4}
\]

So the `x_m=5/48` limit is not an arbitrary successful CFT assignment after the fact.  It is the conformal ground state of the exact periodic dilute-TL module selected by the no-noncontractible-loop condition.

## 4. The central-trinomial state count is the same module dimension

The transparent safe automaton has dimensions

\[
1,3,7,19,51,141,393,1107,3139,\ldots                        \tag{4.1}
\]

for successive widths after the width-one convention is aligned.  These are

\[
T_N=[z^0](1+z+z^{-1})^N.                                     \tag{4.2}
\]

Morin-Duchesne--Kluemper--Pearce's periodic zero-defect standard module has the same central-trinomial dimension sequence.

The earlier note `safe-frontier-dilute-tl-bridge-20260914.md` explains the combinatorial route through vacancies and type-B/no-zero-block annular connectivity.  The present note adds the twist/fugacity information: it is specifically the `alpha=0` member of that zero-defect module which realizes the safe topological semantics.

## 5. Independent descendant check from square-site momentum

The module dictionary makes a sharp prediction for the first descendants of the magnetic primary.  The level-one states

\[
L_{-1}|m\rangle,
\qquad
\bar L_{-1}|m\rangle                                         \tag{5.1}
\]

have

\[
\Delta x=1,
\qquad
s=+1,-1.                                                      \tag{5.2}
\]

The transparent square-site safe transfer was diagonalized independently of the dilute-TL formula.  Its first excited Perron eigenvalue is a twofold degenerate pair.  Restricting the exact one-column translation operator to this eigenspace gives eigenvalues

\[
\boxed{e^{+2\pi i/N},\qquad e^{-2\pi i/N}}                    \tag{5.3}
\]

to numerical precision for `N=5,...,9`.

The scaled energy gap simultaneously obeys

\[
N\log(\lambda_0/|\lambda_1|)\longrightarrow2\pi.            \tag{5.4}
\]

Thus the first relaxation mode is not merely numerically close to `2pi/N`: its lattice momentum is exactly `+/-1`, as predicted by the level-one magnetic descendants of `W_{N,0,i}`.

Machine-readable controls and the reproduction script are in

- `scripts/safe_transfer_momentum_spectrum.py`;
- `results/geometric-consistency/safe-transfer-momentum-spectrum-w5-w9-20260914.json`.

## 6. What this does and does not identify on the square lattice

The square NN/matching safe kernels are not claimed to be Yang--Baxter-integrable dilute `A_2^(2)` transfer matrices.  Their local Bernoulli row weights differ from the integrable triangular-site weights.

The robust identification is instead

```text
annular connectivity module  = dilute periodic zero-defect vocabulary,
noncontractible-loop rule     = alpha=0 / omega=i,
UV continuum ground sector    = magnetic x=5/48,
first descendants             = Delta x=1, momentum +/-1.
```

In other words the square transfer is a nonintegrable lattice regularization acting on the same topological module and flowing to the same `c=0` magnetic sector.

This distinction matters: it licenses the sector dictionary and conformal quantum numbers without pretending that square-site finite-width eigenvalues satisfy the triangular model's Bethe equations.

## 7. Upgrade of the integrable-field-theory programme

`integrable-potts-magnetic-charge-scaling-20260914.md` previously left the critical magnetic twist/source unidentified.  Equation (2.4) removes that ambiguity at the UV lattice/CFT endpoint:

\[
\boxed{\text{target UV sector}=W_{N,0,i},\quad \alpha=0,\quad c_{eff}=-5/4.} \tag{7.1}
\]

The genuinely missing bridge is now narrower:

> Find the **off-critical thermal/massive finite-volume continuation of the `alpha=0`, zero-defect magnetic sector**, or an equivalent excited/twisted Potts NLIE whose UV limit is `W_{0,i}` and whose two thermal signs are related by Potts duality.

This is no longer a search over arbitrary twists.

The 2002 Dorey--Pocklington--Tateo Potts TBA/NLIE gives finite-size thermal flows over continuous `q<=4`, but the specific `q=1` `alpha=0` magnetic-sector finite-volume continuation is not supplied in the material inspected here.  The critical periodic dilute-TL paper fixes exactly which endpoint that continuation must reach.

## 8. Consequence for the dual-odd charge function

The charge scaling function compares two microscopic realizations of the same continuum magnetic module at opposite thermal signs:

\[
\mathcal F(X)=\mathcal E_{W_{0,i}}(X)-\mathcal E_{W_{0,i}}(-X).\tag{8.1}
\]

This makes its leading oddness structurally natural.  It also clarifies the role of the square-lattice sector-odd spin-four correction: the `w^-17/4` mismatch is not a difference of two unrelated primary sectors; it is a lattice-regularization difference **inside the same magnetic module**.

That is precisely why sector-even primary and ordinary irrelevant contributions can cancel so efficiently from the charge root.

## 9. Why higher raw safe eigenvalues should not all be assigned to one Verma tower

The first descendant check is exceptionally clean.  Higher raw eigenvalues of the finite safe Perron block, however, need not form a single irreducible Virasoro Verma module.  Periodic dilute-TL standard modules are reducible/indecomposable at the percolation root of unity, and the finite connectivity space may contain several conformal families and logarithmic extensions.

Therefore the correct next step is **momentum/module resolved spectroscopy**, not fitting every raw energy gap to an integer.

A higher-level identification should be accepted only when its energy, translation momentum and standard-module/character multiplicity all agree.

## 10. Claim boundary

The published periodic standard-module weight formula and the `d=0`, twist-labelled framework are literature facts.  `alpha=0 <-> omega=+/-i` is exact algebra.  The magnetic weight `5/48` follows exactly from the published formula.  The square-site first-descendant momentum/gap check is a deterministic transfer calculation.

For triangular critical site percolation the periodic dilute-TL representation is integrable and explicit.  For square-site NN/matching percolation, the statement is a topological/UV module identification and universality bridge; a full all-width basis/generator intertwiner remains to be proved.
