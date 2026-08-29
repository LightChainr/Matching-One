# Small-width operator spectroscopy: which symmetry actually removes V22?

The first useful transfer-matrix result is a matrix-element zero, not a gap estimate.

We use the exact `Q=4` square-lattice Potts row transfer at its self-dual coupling `v=2`, width four. This is an integer realization of the generic Potts `[2]` carrier; it is not a computation of the percolation threshold or of `x=17/4` from a gap. The 256-state transfer commutes exactly with global `S_4`, row translation, and reflection.

As an unmarked matching insertion, define on a row word the number of monochromatic components on the cycle `C4` minus that on the chord-completed graph `K4`. It is colour blind and invariant under the row dihedral group. Exact central projectors then give, over the full basis,

```text
P_singlet O_matching P_[2,2] = 0.
```

This is the finite-width version of #257's global-endpoint selection rule. The transfer and insertion preserve Potts representation, while the `V_(2,+/-2)` candidate is in `[2]` (realized as `[2,2]` at `Q=4`). Thermal Q4 is a singlet and is not removed.

## What spatial symmetry does—and does not—do

At width four, spin `+/-4` has row-translation phase one. Its reflection-even cosine combination occupies the same row-dihedral sector as a scalar. The joint zero-momentum/reflection-even sector has exact ranks

```text
singlet: 7,
[2,2]:  12.
```

Both sectors are live. Therefore rotation/translation and reflection do not explain the missing lower-dimensional field; at this width they alias it with thermal Q4. Potts charge and the zero-leg/four-leg topology label do the discriminating.

The zero is not an absence theorem. The same colour-blind insertion acts nontrivially inside the charged block:

```text
Tr(P_singlet O_matching) = 4,
Tr(P_[2,2] O_matching)   = 12.
```

So a torus trace or charged two-point function can see `[2,2]`; only the singlet-to-charged endpoint matrix element vanishes. This is the explicit small-width counterexample guarding the claim boundary.

A single colour-transposition seam makes the distinction one-shot. For the unnormalized insertion trace,

```text
singlet twisted/identity = 1,
[2,2] twisted/identity  = 0.
```

Finally, tensoring the insertion with #155's complement-odd marker does not change the Potts projector algebra. Matching parity and Potts colour charge commute: an odd singlet remains unable to couple linearly to `[2,2]`.

The next transfer calculation, if desired, should preserve these explicit projectors and examine the singlet zero-leg block for the thermal Q4 lattice descendant. Increasing width merely to fit a gap would add less information than constructing that operator.
