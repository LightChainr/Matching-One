# Target 1 source-orthogonal continuation

The 2M reveal establishes a bulk `O_far J_D4` response but leaves its complex
transfer phase locked to the JS control. The next question is therefore in
the source plane, not at a larger size and not in qJ.

The marked stream already stores the two Gram entries needed for the unique
least-squares removal of JS from JD:

```text
beta = Re <J_D,J_S> / <|J_S|^2>,
J_D_perp = J_D - beta J_S.
```

Beta is defined separately for each size and orientation at the intrinsic
matching root. It uses no Euler coupling, no transfer and no fitted continuum
exponent. Every leave-one-batch score recomputes both the intrinsic root and
beta after removing that batch. This makes beta a declared source-coordinate
map rather than a coefficient chosen to reduce the result.

The desired coupling is reconstructed linearly from the same rows:

```text
Cov(O_far,J_D_perp)
  = Cov(O_far,J_D) - beta Cov(O_far,J_S).
```

No new acquisition is required. The existing stream does not contain
`|J_D|^2`, so it cannot report the residual source norm or an L2 energy
fraction. That missing statistic does not affect the coupling or its complex
N425/N325 transfer. Orthogonality here refers only to the recorded
same-next-site lattice Gram metric; it is not a claim of CFT field
orthogonality.
