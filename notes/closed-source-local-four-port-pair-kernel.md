# The explicit local four-port pair kernel and its finite Q1 closure

The specified ordered four-port kernel exists and gives the proposed
C4-averaged two-pair indicator exactly. Its single-insertion closure
has a removable Q1 singularity even though its colour-index entries
have poles. It is an additional local invariant interaction direction,
not a change of the old vacant/occupied activity. There is also an
important distinction: C4 averaging preserves orthogonality to the two
old vertex tensors, but does not preserve their full fixed-cut operator
annihilation or the projector property.

Base: `bea717e826df5a22518774b1725ae7bcbe2cb801`. The occupation A has
x vacant throughout this insertion, and its ambient rank, q and E are
kept fixed as requested. No new occupation populations, seam selection,
numerical coupling evaluation, simulation, or source search is used here.

## 1. The rational kernel

At integer Q>=4, let W have orthonormal basis the unordered unequal
colour pairs e={a,b}. Write J_W for the all-ones matrix and
`X_(e,f)=|e intersection f|`. The two lower projectors and their
orthogonal complement are

```
P0 = 2 J_W/[Q(Q-1)],
P1 = [X-4J_W/Q]/(Q-2),
P2 = I-P0-P1
   = I-X/(Q-2)+2J_W/[(Q-1)(Q-2)].                (1)
```

Thus `tr P2=d2(Q)=Q(Q-3)/2`. In particular its entries depend on the
intersection of e and f as follows:

| abs(e intersection f) | P2 entry |
|---:|---:|
| 2 | `(Q-3)/(Q-1)` |
| 1 | `-(Q-3)/[(Q-1)(Q-2)]` |
| 0 | `2/[(Q-1)(Q-2)]` |

Use the specified isometry
`i|{a,b}>=(|a,b>+|b,a>)/sqrt(2)` into ordered colour pairs. Set
`a=N,b=E,c=S,d=W`, with NE incoming and SW outgoing. The kernel
`Pi2=iP2 i^dagger` is

```
Pi2(a,b;c,d) = (1/2)(1-delta_ab)(1-delta_cd) {
  delta_ac delta_bd + delta_ad delta_bc
  -(delta_ac+delta_ad+delta_bc+delta_bd)/(Q-2)
  +2/[(Q-1)(Q-2)] }.                             (2)
```

The factor 1/2 is essential: it comes from the two isometric legs.
This is a globally colour-invariant endomorphism. At integer Q>=4
it is an orthogonal projector of rank d2 on the ordered-pair space.

For a partition pi of the four ports, let D_pi be its equality tensor:
indices in each block are equal, and distinct blocks are not required
to have distinct colours. An expanded partition-diagram formula is

```
Pi2 = (D_NS|EW + D_NW|ES)/2
      -(D_NS + D_NW + D_ES + D_EW)/[2(Q-2)]
      +(D_NES|W+D_NEW|S+D_NSW|E+D_ESW|N)/(Q-2)
      -Q D_NESW/(Q-2)
      +(D_N|E|S|W-D_NE-D_SW+D_NE|SW)/[(Q-1)(Q-2)]. (3)
```

For brevity a single pair such as D_NS leaves the other two ports as
singletons. Formula (3) is an identity of rational equality tensors,
not a rule to recompute topology on its individual virtual joins.

## 2. All 15 external connectivity contractions

Remove the local vacant tensor. Suppose the outside occupied graph
induces connectivity partition pi on the four incident hypergraph
ports. If it has c0 further blocks not meeting these ports, its colour
contraction is

```
Q^c0 F_pi(Q),
F_pi(Q)=sum_(a,b,c,d) Pi2(a,b;c,d) D_pi(a,b,c,d).  (4)
```

The unperturbed local tensor 1 instead gives `Q^(c0+|pi|)`. Here pi is
the exact *graph connectivity*, while D_pi imposes only colour
equalities; coincident colours on disconnected blocks remain summed.

The complete contraction table is

| External partition pi | F_pi | C4-averaged F_pi |
|---|---:|---:|
| N\|E\|S\|W | 0 | 0 |
| NE\|S\|W | 0 | 0 |
| NS\|E\|W | 0 | 0 |
| NW\|E\|S | 0 | 0 |
| ES\|N\|W | 0 | 0 |
| EW\|N\|S | 0 | 0 |
| SW\|N\|E | 0 | 0 |
| NES\|W | 0 | 0 |
| NEW\|S | 0 | 0 |
| NSW\|E | 0 | 0 |
| ESW\|N | 0 | 0 |
| NE\|SW | 0 | d2/2 |
| NS\|EW | d2 | d2 |
| NW\|ES | d2 | d2/2 |
| NESW | 0 | 0 |

This table follows without colour enumeration:

- An input pair NE or output pair SW forced equal is killed by (2).
  Every triple and the all-equal partition contains such a pair.
- The all-free sum is the contraction of P2 with the all-ones vector,
  so is zero. For a single crossing equality, e.g. N=S, it is
  `(1/2) sum_a <v_a,P2 v_a>`, where
  `v_a=sum_(b!=a)|{a,b}>`. The v_a span the incidence/singlet/standard
  subspace eliminated by (1); this contraction is also zero.
- NS|EW closes the ordered-pair identity and gives `tr Pi2=d2`.
  NW|ES closes the swap. The image of Pi2 is symmetric under swap,
  so it gives the same d2.

The algebraic table does not assert that every external partition is
realizable in a given embedded geometry. That is a separate topology
question, not needed for the kernel or its zero contractions.

## 3. C4 average and the actual occupation insertion

Since Pi2 is invariant under a half-turn, the quarter-turn average is

```
Kbar(a,b,c,d) = (1/2)[Pi2(a,b;c,d)+Pi2(b,c;d,a)].   (5)
```

The opposite-pair partition NS|EW is rotation-invariant. The two
adjacent-pair partitions NE|SW and NW|ES are exchanged. Therefore,
relative to the original vacant colour contraction,

```
beta_Q(pi) = (Q-3)/(2Q) {
  I_(NS|EW)(pi) + (1/2)[I_(NE|SW)(pi)+I_(NW|ES)(pi)] }.
                                                        (6)
```

Each I in (6) indicates the exact outside connectivity partition, not
a colour equality test. The ratio uses Q^2 because every surviving
partition has exactly two blocks. All external block factors and the
stipulated original `Q^(-r(A)/2)` cancel; r(A) is not changed by a term
of the diagram expansion.

Replacing the vacant tensor at the named site x by `1+epsilon Kbar`
and leaving the occupied `v delta_all4` term unchanged gives the exact
first-order occupation insertion

```
beta_x(A) = 1_(x vacant) beta_Q(pi_x(A)).           (7)
```

The corresponding numerator inserts the *old* q(A) or E(A) in the
same sum. Unlike the full central seam packet, the local kernel table
alone does not force its support to topology rank1 or its q/E
numerators to vanish. Those conclusions, if true on a specified
geometry, require the outside-connectivity/topology analysis.

## 4. Why the Q1 pole disappears only after closure

For a fixed formal colour pattern the singular part of (2) is

```
Pi2 = -(1-delta_ab)(1-delta_cd)/(Q-1) + O(1).       (8)
```

There is no pointwise finite Q1 kernel obtained by substituting one
colour into those rational entries. In (4), however, the sum of the
numerator of (8) is a polynomial colour count with the factors that
cancel the pole. For the two crossing 2+2 partitions it is `Q(Q-1)`;
for a single crossing equality it is `Q(Q-1)^2`; for all-free ports it
is `Q^2(Q-1)^2`. Other imposed equalities either reduce to these cases
or kill the unequal-pair factor identically. The regular pieces in
(2) give the exact full table, not only its residue.

For example, NS|EW closes directly to

```
Q(Q-1) * (Q-3)/[2(Q-1)] = Q(Q-3)/2.              (9)
```

Thus all 15 contractions are polynomial and have unambiguous finite
values at Q1: 0 or -1 before C4 averaging. The actual ratio (6) is
regular for positive Q, and at Q1 is

```
beta_1(pi) = -I_(NS|EW)
            -(I_(NE|SW)+I_(NW|ES))/2.             (10)
```

The prescription is **contract in the physical outside connectivity
first, then continue its resulting occupation coefficient**. It is
not a literal S1 pair representation. For one insertion this proves
the removable continuation. It does not, by itself, prove that a
uniform finite-strength family with arbitrarily many such insertions
has a pole-free Q1 expansion to every order.

## 5. Which local tangent this is, and which annihilation survives

In the fixed NE-to-SW cut write the original two vertex tensors as
matrices `V0(a,b;c,d)=1` and `V1(a,b;c,d)=delta_(a=b=c=d)`.
For the unaveraged Pi2,

```
Pi2 V0=V0 Pi2=0,
Pi2 V1=V1 Pi2=0.                                (11)
```

The first equality removes the singlet vector, and the second removes
the equal-colour diagonal pair subspace. Their fully contracted
four-tensor overlaps also vanish.

C4 averaging preserves the fully contracted statements

```
sum Kbar V0=0,             sum Kbar V1=0,          (12)
```

as the all-free and all-equal rows of the table show. It **does not**
preserve all of (11) as matrix products in the original cut. For
distinct colours a,b,

```
Kbar(a,a,b,b) = (Q-3)/[4(Q-1)] != 0,   Q>=4.      (13)
```

Thus Kbar acting on an equal-colour input |aa> has a nonzero equal-
colour output |bb>; it does not annihilate V1 as an operator. The
rotated summand has recoupled which two ports form the pair. Kbar is
a C4-invariant tensor, not the same fixed-cut orthogonal projector.
Indeed the normalized vector `(|aa>-|bb>)/sqrt(2)` has quadratic form
`-(Q-3)/[4(Q-1)]<0` under the real symmetric Kbar, so Kbar cannot be
an orthogonal projector. The quarter-turn here reshuffles input and
output ports; it is not unitary conjugation within the old ordered-pair
space. This does not invalidate (6), (7), or the local perturbation.

The old thermal/activity family lies in `span{V0,V1}`. At integer
Q>=4, (12) and the nonzero table prove that Kbar is outside that span
under the ordinary positive Frobenius inner product. Accordingly:

- It is a genuine additional **local colour-vertex interaction tangent**
  of `V0+epsilon Kbar+v V1`. For fixed integer Q>=4, sufficiently small
  epsilon retains nonnegative vacant tensor entries, since V0 is
  strictly positive and Kbar has finitely many finite entries.
- It is not a thermal reparametrization, a normalization of V0/vV1,
  or a tangent already generated by changing the old occupancy p.
- Its Q1 continuation is the bounded geometric insertion (7),(10),
  not an entrywise positive pair-colour tensor at Q1. A nonthermal
  microscopic direction does not alone prove nonzero response of U.

Finally, the endpoint identity `ell P_[2]=0` concerns an invariant
single endpoint into the trivial representation. Here the two legs
of W and W* contract through an endomorphism; the identity/swap
closures have value d2. The insertion therefore does not contradict
that endpoint rule. It is also not the full central **seam** projector:
it acts on four ports at a named vacant vertex and selects the local
outside 2+2 connectivity. The global seam acts on all winding colour
lines with a class sum and has the previously proved rank1 support.
Equating these two contractions would require an additional identity,
which the present kernel calculation neither assumes nor supplies.

## Scientific card

The new mechanism is an explicit local invariant vertex perturbation
whose closed Q1 statistic is a fixed two-pair connectivity indicator.
All Bell4 closures and their coefficients are determined. A removable
entrywise pole leaves a finite signed local response; the nontrivial
trace does not pass through the forbidden invariant single endpoint.
The C4 recoupling caveat (13) prevents mistaking this local interaction
for a pure projector propagation in one fixed pair channel. No claim
about a continuum four-leg field, Jordan structure, or a nonzero
original-U response follows without the separate topology and
transmission calculations.
