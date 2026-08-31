# The finite [22] trace enters original U through its normalizer

## One extra positive seam isolates a definite colour sector

At the already specified m=2, Q=4 local colour law, keep the **same spatial
identity colour seam**, geometry, activity `a=y/32` and local weights in
all three partitions. Let

```text
T = Z(1,1),
D2 = Z(1,(12)(34)),
C3 = Z(1,(123)).
```

Each is a nonnegative local partition function; the empty configuration
makes each strictly positive at finite positive activity. The spatial
identity commutes with every temporal colour permutation. Any geometric
row-translation seam is kept fixed and commutes with these colour actions.

The exact S4 central-projector formula gives

```text
Z22 = Tr(P_[22] transfer) = (T+3 D2-4 C3)/6.                 (1)
```

The factor is `dim([22])/|S4|=2/24`; class sizes are1,6,3,8,6 and
the [22] characters are2,0,2,-1,0. In particular (1) eliminates every
other S4 irrep, not just a chosen singlet competitor. The complete
relative temporal-character table is

| irrep | dimension | transposition | double transposition | 3-cycle | 4-cycle |
|---|---:|---:|---:|---:|---:|
| [4] | 1 | 1 | 1 | 1 | 1 |
| [31] | 3 | 1/3 | -1/3 | 0 | -1/3 |
| [22] | 2 | 0 | 1 | -1/2 | 0 |
| [211] | 3 | -1/3 | -1/3 | 0 | 1/3 |
| [1111] | 1 | -1 | 1 | 1 | -1 |

The [existing five-partition reconstruction](https://github.com/LightChainr/Matching-One/blob/977fea9272c780aea19cc47f8d33324c28a1293e/notes/closed-source-hypergraph-rc-twist-projection.md)
already contains T and D2. C3 is one additional partition **type**, not
a claim that the existing five numerical p-jets have been acquired.
It supplies information beyond the three total rank weights.

## Why the original five seams are not already a colour-sector test

At m=2 their colour translations form the Klein four group V4. The [22]
representation is trivial on V4, as are both one-dimensional irreps.
In the spatially untwisted sector the identity/double-transposition traces
therefore cannot separate their multiplicities. Formula (1)'s 3-cycle can.

The five partitions also do not all use the same spatial sector. Writing
`s=(12)(34), t=(13)(24)`, they separate into

```text
spatial1: Z(1,1), Z(1,s);
spatial s: Z(s,1), Z(s,s), Z(s,t),
```

where the last partition is colour-conjugate to the old full-rank twist.
The second row has centralizer D8, not full S4, as its boundary symmetry.
Although a temporal transposition commuting with s is legal and positive,
its trace cannot simply be assigned the untwisted S4 character ratio0.
Indeed `[22] restricted to D8 = trivial + restricted sign`; the two
components can have different transfer weights. Mixing those spatial
sectors into one character test would give a false discriminator.

## The central projection is supported only on rank1

For an occupation configuration A put `w(A)=a^K 4^(B+C_B)` in the
untwisted colour law. Its occupied NN homology is saturated, as in the
cited twist construction. At rank0 all three seams have the same colour
count, so (1)'s contribution is zero. At rank2 there is exactly one
rank2 essential component: two disjoint essential components with
independent slopes would intersect on the torus. Its homology is Z²;
all other components are contractible. The double transposition has no
fixed colour and the 3-cycle has one, hence `D2(A)=0`, `C3(A)=w(A)/4`.
Again the contribution to (1) is zero.

At rank1 let c be the number of essential occupied components. Their
common primitive deck slope is `(u,n)`, where n is the **temporal deck
winding**, not Cartesian displacement. A component winding n times sees
the permutation to the nth power. Thus exactly

```text
Z22(A) = w(A) f22(n,c),
f22(n,c) = {1+3 1_(2 divides n)
             -4[1_(3 divides n)+1_(3 does not divide n) 4^(-c)]}/6. (2)
```

Set f22=0 at ranks0,2. For example f22=0 for n=1,c=1, but f22=1/8
for n=1,c=2. For n=3 mod6 it is -1/2. This **signed character component
is not a probability or a nonnegative population share**. Neither rank
alone nor the old `(K,g,q)` histogram determines it.

This support theorem does not assert that the trace component vanishes.
On an axis LxL torus, L>=5, its activity polynomial begins

```text
Z22(a)=L(L-3) 4^(2L) a^(2L)+O(a^(2L+1)).                    (3)
```

To attain the first term, select two nonadjacent vertical occupied columns.
There are L(L-3)/2 choices, c=2, n=1, B=K=2L, and each contributes
`2*4^(2L) a^(2L)` to (1). Below2L sites there cannot be two essential
components. A single component with nonzero even n requires a primitive
slope with length at least3L; an odd multiple of3 costs at least4L.
A single n=±1 component contributes zero; n=0 also contributes zero.
At2L sites the only contributing two-component configurations are those
two straight vertical columns. This proves (3), including the absence
of cancelling lower terms, without enumerating configurations. For N25
the leading integer coefficient is `10*4^10=10485760`.

## Exact normalizer-mediated interface to the original q/E and U

Let L_j be the untwisted occupation partition restricted to rank j.
Use two explicit normalizations:

```text
Wstar = L0 + L1/2 + L2/4,
mathcal D = T+R = 4 Wstar,
Wstar22 = Z22/2,
mathcal D22 = 2 Z22,
z22 = Wstar22/Wstar = 2 Z22/mathcal D.                       (4)
```

All configuration-independent prefactors have been removed. The last
two factors must not be interchanged. Because q=r-1 and E=q² both vanish
at rank1, their **unnormalized [22] numerators are exactly zero**.

Multiply this fixed central component by a bookkeeping parameter x,
`Wstar(x)=Wstar+(x-1)Wstar22`, while its zero q/E numerators stay zero.
For each geometry, at fixed p,

```text
q_g(x)=q_g/[1+(x-1)z22_g],
E_g(x)=E_g/[1+(x-1)z22_g],
j_x q_g=-q_g z22_g,       j_x E_g=-E_g z22_g.                 (5)
```

The character component is signed; x labels this fixed algebraic
attribution, not a new fitted colour number or a claimed local positive
source. Equation (5) is well-defined near x=1, where Wstar is positive.

Now retain the actual separate geometry normalization, pooled root and
thermal slope. Set `M=mean(q_g)`, `Y=P4(E_g)`, `D=M_p`,
`U=A_N Y_p/D`, `hM=-mean(q_g z22_g)`, `hY=-P4(E_g z22_g)`.
At the original pooled M=0 root,

```text
J22 = d_x U |x=1
    = A_N [ hY_p/D - Y_pp hM/D²
             -Y_p hM_p/D² + Y_p M_pp hM/D³ ].               (6)
```

This supplies the finite-colour matrix-element interface that the
[regular endpoint zero](weak-q-paths-and-regular-selection.md) alone
could not determine: the raw q/E numerator is blind, but the normalizer
has the concrete nonzero trace polynomial (3), and its transport to U
is (6). There is no identity here forcing J22=0. **Its value at the
specified matching root has not been computed**, and (3) alone does not
establish a nonzero J22.

At a single geometry's q=0 root the denominator change cannot move that
root. At a two-geometry pooled root q1+q2=0 need not mean q1=q2=0;
if z22_1 differs from z22_2, `p_x=mean(q_g z22_g)/D` can be nonzero.
Dropping this term would change the original U interface.

## The next finite decision is now specified

For the same N25 pair and m=2, obtain the fixed central component's p-jet
and evaluate (6) once. A strict nonzero bound would reject complete
normalizer-neutrality of this [22] component for that finite U; zero
would not erase its trace polynomial or establish continuum neutrality.
Existing `(K,g,q)` counts supply the baseline but lack the `(n mod6,c)`
rank1 information in (2). Obtain that specified information or compute
the three compatible seam partitions, without fitting a new descriptor.

The measurement definition is [machine-readable](../analysis/four_leg_trace_interface.json).
No new data, transfer diagonalization, enumeration or test campaign was
performed for this derivation. It identifies a finite Q=4 colour channel,
not uniquely the continuum V_(2,2) field. It is not a Q->1 derivative or
a confluent continuation; those require their own explicit continuation.
The independent P154/P334/F4 decisions and the fixed weak-Q path control
remain unchanged.
