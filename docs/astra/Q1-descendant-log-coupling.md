# Q1 — Does the `h = 5/8` logarithmic coupling descend to the level‑4 `x = 21/4, s = 4` state?

**Blocks:** issue #275, currently the project's only P0, closed as
`UNIDENTIFIABLE_WITH_CURRENT_ASSETS` for a reason that no further sampling can fix.

---

## The question

In the `c → 0` bulk CFT of two-dimensional critical percolation, the thermal
(energy) sector at `h = 5/8` carries a rank‑2 Jordan pair whose logarithmic
coupling is known and equal to `μ = -5/4` in the normalization of
Y. He, *Logarithmic operators in c=0 bulk CFTs*, arXiv:2411.18696 (SciPost Phys. 19, 008).

We use a construction that carries this pair to a **level‑4 generalized state with
`x = 21/4` and spin `s = 4`**, through a non-null thermal `Q4` descendant. What
survives the lift is the affine scale cocycle: writing `D` for the dilation
generator, the image pair still satisfies

```text
D q = x q,        D q~ = x q~ + q,
```

so the rank‑2 structure and its `log L` law propagate.

> **Is the logarithmic coupling of that level‑4 pair determined by the primary-level
> normalization — and if so, what is its value?**

Equivalently: is the ratio fixing the two-point normalization of `(q, q~)` at
`x = 21/4, s = 4` a function of `μ = -5/4` and the descent data alone, or is it a
free parameter of the construction?

## Why this is not automatic

For an ordinary Virasoro descendant the two-point normalization is forced by the
algebra, and the question would be a computation. Here the descent is through a
`Q4` structure of the `Q`-state Potts / Temperley–Lieb setting rather than a plain
Virasoro lowering, and the interchiral structure at `c → 0` is exactly where the
usual argument can fail. Our own note records that we have **not** established

- that the lattice matching involution is the required interchiral automorphism,
- that this particular level‑4 generalized state has nonzero lattice coupling, or
- that its logarithmic coupling is fixed by the primary-level normalization.

The third is this question. The first two are separate and are not being asked here.

## Why the answer decides something

The observable we are trying to identify is scored by a profile-rank calculation
against a frozen covariance. Two candidate mechanisms are on the table: an
**ordinary semisimple** correction with a transfer parameter `κ`, and the
**Jordan / logarithmic** image above. The design calculation gives

```text
rank(semisimple, κ fixed) = 8
rank(Jordan)              = 8
rank(combined)            = 12
dim(image intersection)   =  4
```

so the two are not separated: **with unrestricted amplitudes, the closure of the
semisimple family as `κ → 1⁻` contains the whole Jordan image.** That is the
familiar collision — a semisimple pair degenerating into a Jordan block — and it is
why the verdict is `PARTIALLY_IDENTIFIABLE`, not identification.

The degeneracy is a statement about *unrestricted amplitudes*. If the logarithmic
coupling of the level‑4 pair is fixed, the Jordan candidate's amplitude is no longer
free, its image is a proper subvariety rather than a subspace, and the `κ → 1`
closure need no longer contain it.

## Decision rule

| Answer | What we do |
|---|---|
| **Fixed, with a value** | The value becomes the candidate-specific column. We score it once against the existing frozen covariance (no new data), and #275 reopens with a real test. |
| **Fixed in form but not in value** (determined up to one unknown constant) | We ask which *ratios* across two moduli are independent of that constant, and score those. This is still a reopening. |
| **Genuinely free** | The `κ → 1` degeneracy is structural, not an artifact of our design. `UNIDENTIFIABLE_WITH_CURRENT_ASSETS` is upgraded from a status to a conclusion, the Q4/Jordan identification line is closed rather than parked, and effort moves to a different observable. |

The third outcome is worth the query. Closing a line on a proof is a result; leaving
it parked indefinitely is not.

## Context you may assume, and need not re-derive

- Percolation as `Q → 1` of the `Q`-state Potts model; `c = 0`; the standard
  Kac table and the `h = 5/8` thermal identification.
- The `c → 0` catastrophe and the general theory of rank‑2 Jordan cells and their
  `b`-numbers, including the stress-tensor case.
- `μ = -5/4` at `h = 5/8` per arXiv:2411.18696, and the companion
  arXiv:2510.09868.
- Jacobsen–Ribault–Saleur, arXiv:2208.14298, for the torus sector decomposition
  conventions.

## Do not spend output on

- Re-deriving `μ = -5/4`, or reviewing `c → 0` logarithmic CFT in general.
- The physical interpretation of `x = 21/4` or of spin‑4 observables.
- Whether the lattice model actually couples to this state — that is a different,
  known-open question and we are not asking it.
- Numerical estimates of anything. This question has no numerical content.
- Advice about experiment design, sample sizes, or statistics.

A three-line answer with a reason is worth more here than a survey.

## Provenance of the framing above

The construction and the three unestablished points are quoted from this
repository's own notes (`notes/thermal-q4-jordan-inheritance.md`,
`notes/thermal-jordan-spin4-descendant.md` §7). The rank and intersection figures
are from the frozen forward-identifiability score recorded on issue #275. If any
part of the framing is wrong, say so — a corrected framing is a useful answer.
