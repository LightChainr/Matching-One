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

## 5. Minimal toroidal rank-source trace dictionary

The same-parameter topological count algebra elsewhere in this project has the schematic form

```text
G(s,t)
 = endpoint_rank2 * s
   + endpoint_rank0 * t
   + neutral_rank1 * H(st).                              (5.1)
```

The loop dictionary shows exactly how the neutral product variable should enter a toroidal modified trace:

```text
z=st = alpha^2/Q,
alpha=sqrt(Q s t).                                        (5.2)
```

The two endpoint variables cannot be represented by `alpha`, because both rank0 and rank2 sit in the zero-noncontractible-loop sector. Introduce two zero-through-line amplitudes, schematically

```text
q0 : trivial/rank0 topology,
q2 : cross/rank2 topology.
```

Then the minimal toroidal source object has the typed structure

```text
T(q0,q2,alpha)
 = q0 Z0
   + sum_u Z_{1,u}(alpha^2/Q)
   + q2 Z2.                                               (5.3)
```

After slope aggregation and normalization, the percolation rank-source specialization is

```text
q0=t,
q2=s,
alpha^2=st                 (Q=1).                         (5.4)
```

Thus the exact finite rank-source algebra itself predicts the minimal number and roles of continuum seam/trace parameters:

- one continuous noncontractible-loop fugacity for the neutral gas;
- two endpoint amplitudes, or equivalently one normalization plus one odd endpoint/duality source.

This is a target interface, not a claim that the massive theory has already supplied these three parameters.

## 6. Stronger #782 target

This reduces the continuum problem further. A useful massive affine-TL object would be a toroidal modified trace

```text
T_massive(mL,mR; q0,q2,alpha)
```

such that:

1. `alpha` weights noncontractible loops and hence evaluates `H(z)` with `z=alpha^2/Q` in rank-one sectors;
2. `q0,q2` separate the two zero-through-line topologies (rank0 and rank2);
3. the UV limit reproduces the Arguin/Pinson critical homology weights and the finite rank-source algebra;
4. the IR one-particle/loop terms reproduce the correct complete winding insertion, not merely an open-boundary crossing amplitude.

This is more specific than asking for a generic twisted Potts TBA.

## 7. A finite-width test before any continuum construction

On an existing periodic FK/TL transfer matrix, introduce a symbolic or numerical noncontractible-loop fugacity `alpha` while leaving contractible loops at `sqrt(Q)`. In a rank-one sector, compare the transfer output to a direct essential-component count generating function under

```text
z=alpha^2/Q.
```

The coefficients must agree configuration by configuration after the standard toroidal FK/loop normalization is included. A failure would indicate a convention/Markov-trace mismatch before any massive interpretation is attempted.

A second finite test is to turn on independent zero-through-line endpoint amplitudes and verify (5.3) against directly classified rank0/rank1/rank2 configurations.

## 8. Claim boundary

- Exact within a fixed rank-one FK homology sector: `N_nc=2K` and `z=alpha^2/Q`.
- Exact source typing: the neutral variable belongs to `alpha`; trivial/cross separation requires an additional zero-through-line projector.
- Structural conclusion: a noncontractible-loop seam naturally realizes the neutral-count source and is better suited to Q-continuation than literal spin permutations.
- Not claimed: an existing massive formula for the zero-through-line projector, a completed Q->1 TBA, or direct transfer of this bond-FK loop identity to square-site complete-component amplitudes without a universality/sewing argument.
