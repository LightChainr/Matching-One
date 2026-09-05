# Q3 — What is the smallest escape from the Gaussian-unit no-go?

**Blocks:** the Gaussian-cover production line, which a theorem of ours has just
declared incapable of answering the question it was built for.

---

## The no-go we proved

We work with percolation on planar quotients `Z[i]/(g)` — square-lattice tori
realized as Gaussian ideal quotients, chosen because they carry exact deck-character
and cover arithmetic (norm 2 gives `Z₂`, norm 5 gives `Z₅`, norm 4 via `2i` gives
`Z₂ × Z₂`, norm 10 gives `Z₁₀`, and `(1+i)² = 2i` gives an exact coarse/detail
Hadamard split).

The observable is a `C₃` homology character: three probabilities `P₀, P₁, P₂` of
distinct unoriented homology target lines, combined as

```text
z = (P₀ - π₀) + ω²(P₁ - π₁) + ω(P₂ - π₂),        ω = e^{2πi/3},
```

against a continuum baseline `π_j`. The intent was to read a local spin from the
phase of `z`, distinguishing spin‑0, spin‑4 and spin‑8 transport.

**The theorem.** For every Gaussian ideal quotient, multiplication by the unit `i` is
an automorphism of both the graph and the bond measure. In the period basis `(g, ig)`
it acts on homology as `(m, n) ↦ (-n, m)`, which exchanges the unoriented target
lines `l₀` and `l₁`. Hence

```text
P_{l₀} = P_{l₁}     at every finite norm and every bond probability p,
```

and the `τ = i` continuum baseline obeys the same equality. Therefore

```text
z = ω · ((P₂ - π₂) - (P₀ - π₀))  ∈  ω·R.
```

The phase is pinned to a real line for the entire family, exactly and non-asymptotically.
Signed-real spin‑0 transport is a *theorem* here, not a measurement. Our earlier
spin‑8 reading at one geometry was an angle alias caused by omitting the spin‑0
hypothesis at an angle where `-8δ` happened to be within a degree of `π`.

Consequence: **no further production in this family, with this observer, can
identify local spin.** One must leave the unit-rotation symmetry class or change the
observer.

## The question

> Characterise the planar lattice/quotient families with exact deck-character
> arithmetic that admit an observer whose `C_n` homology character is **not** forced
> onto a real line by a unit (or other automorphism) of the underlying ring — and
> give the smallest concrete example.

Three directions we can see, and would like judged rather than merely listed:

1. **Change the ring.** Eisenstein integers `Z[ω]` have unit group of order 6, which
   looks worse, not better — but the triangular geometry pairs naturally with a `C₃`
   character, so the interaction of the two is not obvious to us. Is there a ring of
   integers, or a class of non-principal ideals, where the unit action does *not*
   act transitively on the relevant homology lines?
2. **Change the observer.** Keep `Z[i]` and break the `i`-symmetry in the observable
   rather than the lattice — a decorated, marked or sublattice-restricted target
   whose homology classes are not permuted by `(m,n) ↦ (-n,m)`. Is there such an
   observer that still has a computable continuum baseline?
3. **Break the measure.** An anisotropic bond measure destroys the `i`-automorphism
   immediately, but also destroys the isotropy that makes the continuum baseline
   meaningful. Is there a controlled anisotropy — a marginal direction — that breaks
   the automorphism without moving the fixed point?

If a general obstruction rules out all three, that is the most valuable answer.

## Decision rule

| Answer | What we do |
|---|---|
| **A concrete family/observer escapes** | Build it. The exact cover arithmetic is the reason this line existed; if it can be kept while breaking the unit symmetry, the whole production programme restarts on the new family. |
| **Escape exists but loses exact cover arithmetic** | We weigh a computable baseline against exact deck characters, which is a real trade and one we can only make once someone tells us it is the trade. |
| **General obstruction** | The line closes on a theorem rather than on fatigue, and the exact results already obtained stand as a negative structural result rather than as a stalled programme. |

## Context you may assume

- Percolation on tori, homology classes of crossing clusters, and the standard
  torus-sector conventions of Jacobsen–Ribault–Saleur, arXiv:2208.14298.
- Arithmetic of `Z[i]`, `Z[ω]`, ideal quotients, deck groups of abelian covers.
- That we can compute exactly on finite quotients — configuration enumeration at
  small sizes and Monte Carlo at large ones. Feasibility of computation is not the
  constraint; the symmetry is.

## Do not spend output on

- Re-proving the no-go above. We hold it as exact and it has been checked by
  direct enumeration.
- General exposition of percolation on the torus.
- Estimating how large a system would be needed for anything.
- Suggesting that we simply use a different observable without saying which — the
  content of this question is entirely in the specific construction.

## Provenance of the framing above

The theorem is our commit `74f55006`, recorded on issue #275; the cover-arithmetic
list is `docs/RESEARCH-MAP.md` §D. The angle-alias history is on the same issue and
is included only to explain why the observer, not the statistics, is the problem.
