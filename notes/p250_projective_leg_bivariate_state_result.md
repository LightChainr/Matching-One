# P250/P249/P255 result: C4 survives, low-dimensional common state does not

The fresh N505 mixed-displacement stream gives a clean separation between
symmetry and state dimension.  C4 covariance is present; a shared commuting
diagonal state of rank one, two, or three is not.

## Exact structural result before statistics

The existing `(parent,fiber)` section carries a position-dependent C4 cocycle:

```text
j' = k j + s(x) mod 5,
k_plus=3, k_minus=2.
```

The exact gate constructs a gauge satisfying
`t(Rx)=s(x)+k t(x) mod 5` on all 101 parent sites, with zero failures.  The
fresh stream therefore measures a genuinely C4-covariant charged pair rather
than asking a constant matrix to absorb a section artifact.

The signed C4 orbit score then closes: `chi2=66.858/64`, `p=0.379`.  The
rank-two and rank-three x/y characteristic polynomials are also compatible
(`p=0.851` and `0.960`).  Ordinary directional symmetry breaking is not the
reason the prior scalar model failed.

## Frozen rank result

Only the two axes through degree three fit each candidate rank.  The first
mixed points `xy,x^2y,xy^2` test the common commuting product; all total-degree
four points are held out.

| rank | mixed gate | degree-4 heldout |
|---:|---:|---:|
| 1 | `1203.36/24`, `p=4.72e-239` | `335.67/40`, `p=2.23e-48` |
| 2 | `221.67/24`, `p=6.32e-34` | `150.86/40`, `p=8.98e-15` |
| 3 | `66.98/24`, `p=6.21e-6` | `141.73/40`, `p=2.69e-13` |

Thus the minimal rank is not two and not three: the frozen decision is
`no_commuting_common_rank_le_3`.  The data have the correct C4 algebraic
envelope but require a richer internal realization to reproduce mixed
translations.

This does **not** say exact lattice translations fail to commute.  The rejected
model is a low-dimensional simultaneous diagonal exponential realization.
Live explanations include rank at least four, Jordan/non-diagonal structure,
finite-torus images, or a context-enriched projected state.

## Scientific card

- **Mechanism space changed:** the two-state interpretation is removed; even
  rank three cannot transport the axis fit into mixed displacements.
- **Not proved:** a unique higher rank, noncommutation, Jordan form, or a
  continuum field count.
- **Observer/sector/source/geometry:** C4-gauged neutral projective-leg pairs;
  charges 1/2, plus/minus Gaussian children; fresh N505 radius-four diamond.
- **Dependency group:** new seed `25050510120261130`, counters `[0,80000)`;
  independent of the earlier N505 axis stream.
- **Next discriminator:** use this same diamond first for a frozen,
  covariance-whitened block-Hankel rank lower bound.  Only if transfer matrices
  of rank four to six remain needed should the grid extend to degree five for
  an independent shifted holdout.

`D^5` is not tested: simultaneous deck phases cancel exactly in this neutral
two-point row.
