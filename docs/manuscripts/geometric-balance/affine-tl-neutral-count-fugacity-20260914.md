# Affine Temperley--Lieb fugacity for the rank-one neutral count

Date: 2026-09-14

Status: exact topological/loop-weight dictionary in a fixed rank-one FK homology sector, plus a narrowed massive-continuum interface for #782.

## 1. Rank-one FK topology and medial noncontractible loops

Consider a toroidal FK configuration in a fixed primitive rank-one homology sector. Let

```text
K = number of essential FK clusters,
```

all of which are parallel to the same primitive homology line.

Take disjoint regular neighbourhoods of the essential clusters. Each rank-one essential cluster is an annulus with contractible holes attached. Its two annular boundary components are noncontractible and parallel to the cluster homology. Contractible holes add only contractible boundary loops.

Because distinct essential FK clusters are disjoint and parallel, their annular boundary components are distinct. Therefore the number of noncontractible medial/interface loops is exactly

```text
boxed: N_nc = 2 K.                                       (1.1)
```

This is a topological statement within the rank-one sector; branch decorations and contractible holes do not change `N_nc`.

## 2. Change the noncontractible loop fugacity

In the standard FK-to-loop representation, the ordinary loop fugacity is

```text
n = sqrt(Q).
```

Now keep the contractible-loop weight equal to `sqrt(Q)` but assign a separate weight

```text
alpha
```

to each noncontractible medial loop. Relative to the ordinary FK weight, a rank-one configuration with K essential clusters is multiplied by

```text
(alpha/sqrt(Q))^(N_nc)
 = (alpha^2/Q)^K.                                        (2.1)
```

Hence, for the slope-resolved neutral-count generating function

```text
Z_{1,u}(z)
 = sum_{rank1 slope u configs} ordinary_FK_weight * z^K,
```

the affine/periodic TL noncontractible-loop fugacity realizes it exactly through

```text
boxed: z = alpha^2/Q.                                    (2.2)
```

This dictionary is independent of integer Potts spin permutations and therefore has a natural analytic continuation in Q, at least at the algebraic loop-model level.

At percolation `Q=1`, simply

```text
z=alpha^2.                                                (2.3)
```

So the common-window neutral law `H(z)` is not merely analogous to an affine-TL seam response: in a fixed rank-one sector, its fugacity is precisely the square of the noncontractible-loop weight at Q=1.

## 3. Relation to the permutation-twist dictionary

For integer Q, `potts-permutation-twist-homology-20260914.md` gives

```text
z_u = |Fix(sigma^a tau^b)|/Q
```

for a rank-one slope `u=(a,b)`.

Combining with (2.2), the corresponding effective loop fugacity is

```text
alpha_u = sqrt( |Fix(sigma^a tau^b)| ).                   (3.1)
```

Thus the permutation-twist evaluation points are a discrete subset of the continuous affine-TL fugacity line. The loop/seam variable supplies the analytic interpolation that literal `S_Q` permutations cannot provide as `Q->1`.

## 4. What alpha does not resolve

The same topological reason that makes (2.2) simple also exposes the remaining blocker.

For rank 0 and rank 2/cross topology, the standard medial-loop classification has no family of parallel noncontractible boundary loops whose count is `2K`. Changing only `alpha` therefore primarily resolves the rank-one neutral gas; it does **not** by itself distinguish the trivial and cross endpoint sectors.

That distinction is carried by the toroidal Euler/duality correction in the random-cluster/loop correspondence and by the special zero-through-line amplitudes of the modified Markov trace.

Therefore a complete massive rank-source construction needs two typed ingredients:

```text
alpha (or z=alpha^2/Q) : neutral rank-one count/source,
q_top                    : zero-through-line trivial-vs-cross projector.   (4.1)
```

The second symbol is schematic here; no existing massive formula is asserted.

## 5. Stronger #782 target

This reduces the continuum problem further. A useful massive affine-TL object would be a toroidal modified trace

```text
T_massive(mL,mR; alpha, q_top)
```

such that:

1. `alpha` weights noncontractible loops and hence evaluates `H(z)` with `z=alpha^2/Q` in rank-one sectors;
2. `q_top` separates the two zero-through-line topologies (rank0 and rank2);
3. the UV limit reproduces the Arguin/Pinson critical homology weights;
4. the IR one-particle/loop terms reproduce the correct complete winding insertion, not merely an open-boundary crossing amplitude.

This is more specific than asking for a generic twisted Potts TBA.

## 6. A finite-width test before any continuum construction

On an existing periodic FK/TL transfer matrix, introduce a symbolic or numerical noncontractible-loop fugacity `alpha` while leaving contractible loops at `sqrt(Q)`. In a rank-one sector, compare the transfer output to a direct essential-component count generating function under

```text
z=alpha^2/Q.
```

The coefficients must agree configuration by configuration after the standard toroidal FK/loop normalization is included. A failure would indicate a convention/Markov-trace mismatch before any massive interpretation is attempted.

## 7. Claim boundary

- Exact within a fixed rank-one FK homology sector: `N_nc=2K` and `z=alpha^2/Q`.
- Structural conclusion: a noncontractible-loop seam naturally realizes the neutral-count source and is better suited to Q-continuation than literal spin permutations.
- Not claimed: an existing massive formula for the zero-through-line projector, a completed Q->1 TBA, or direct transfer of this bond-FK loop identity to square-site complete-component amplitudes without a universality/sewing argument.
