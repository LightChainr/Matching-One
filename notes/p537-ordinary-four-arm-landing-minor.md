# An ordinary four-arm landing block already has a nonzero Schur minor

## Decision

The proposed rank-one cancellation theorem is false at finite landing
level.  On the ordinary `4 x 4` NN torus there are three physical fibres,
each with exactly two opposite occupied and two opposite vacant incident
arms at the thermal site, for which the C4-averaged, root-conditioned
source/readout block has the exact minor

\[
                              \boxed{-\chi/2\ne0}.             \tag{1}
\]

Here `chi` is the nonzero angular coefficient multiplying `E` in the
chosen P4 row.  The result is independent of the root coefficient `R`, of
the Schur coefficient `beta`, and of the common Bernoulli centering.

Thus an additional fifth or sixth arm is **not** structurally necessary
for a nonzero projected landing channel.  This stops the route

```text
ordinary four-arm block is rank one
  => three-packet term is only a thermal reparametrization
  => first surviving term needs four packets.
```

It does not prove that the resulting infinite-volume signed spatial sum is
nonzero or large.  That remains a probability/transport problem, not a
finite-algebra identity.

## 1. The projected local coordinates

For a fixed pair `(x,y)` and thermal site `z`, retain the exact audit
coordinates

\[
 k=D_z g_{xy},\qquad
 h_R=D_z(y-Rm)=\chi D_zE-RD_zq.                  \tag{2}
\]

The first coordinate is the kernel-reconnection landing channel.  The
second is the readout pivotal after quotienting by the matching/root
thermal direction.  Entry and completion have

\[
\begin{array}{c|cc}
 &D_zq&D_zE\\ \hline
 0\to1&1&-1\\
 1\to2&1&+1,
\end{array}
\quad\Longrightarrow\quad
 h_{entry}=-\chi-R,\qquad h_{completion}=\chi-R. \tag{3}
\]

A rank-preserving reconnection has `D_zq=D_zE=0` but may have `k!=0`.
These are the two exact terms in the midpoint covariance derivative; they
are not two new sources.

The remaining Schur term `-beta B` depends only on the Bernoulli thermal
score at this finite fibre.  We choose all three backgrounds with the same
`N=16` and `K=6`, so it adds the same scalar to their source coordinate.
Taking landing differences removes it exactly.  No value of `beta` is
assumed.

## 2. Three physical ordinary-four-arm fibres

Coordinates below are modulo four.  All unlisted sites are vacant, as are
`x,y,z`.  At every `z`, its incident `(N,E,S,W)` occupation is either
`(1,0,1,0)` or `(0,1,0,1)`.  Hence there are exactly four alternating
incident arms and no local fifth branch.

### Entry fibre `A`

```text
occupied = {(0,0),(1,0),(0,1),(2,0),(1,1),(0,2)}
z=(3,0),       (x,y)=((2,1),(1,2))
rank: 0 -> 1,  16 g_xy: 0 -> 0.
```

Filling `z` closes the almost-complete horizontal carrier.  Thus

\[
                         (k,h_R)_A=(0,-\chi-R).               \tag{4}
\]

### Kernel-only fibre `B`

```text
occupied = {(0,0),(1,0),(0,1),(2,0),(1,1),(2,2)}
z=(2,3),       (x,y)=((2,1),(1,2))
rank: 0 -> 0,  16 g_xy: 4 -> 0.
```

The opposite occupied arms merge two zero-image components and reroute the
two marked four-port signatures without changing ambient rank.  Therefore

\[
                         (k,h_R)_B=(-1/4,0).                  \tag{5}
\]

### Completion fibre `C`

```text
occupied = {(0,0),(1,0),(0,1),(2,0),(0,2),(3,0)}
z=(0,3),       (x,y)=((1,1),(2,1))
rank: 1 -> 2,  16 g_xy: 4 -> 4.
```

The horizontal carrier is already essential; filling `z` closes the
independent vertical carrier.  Hence

\[
                         (k,h_R)_C=(0,\chi-R).                 \tag{6}
\]

The three states distinguish the actual landing identification of the two
occupied arms: closing the first essential cycle, merging two contractible
branches, or closing the second cycle.  No arm-number statistic can identify
these three outcomes.

## 3. C4 and Schur leave an exact nonzero minor

Rotate each complete tuple `(A,x,y,z)` through its C4 orbit.  Ambient rank,
the canonical `Kreg` value and the entry/completion type are unchanged.
The spin-four character is `+1` under a quarter turn, so orbit averaging
leaves (4)--(6) unchanged rather than cancelling them.

Subtract fibre `A` from `B` and `C`.  Common source centering and the
common `-beta B` Schur column disappear.  The resulting source/readout
matrix is

\[
 \mathcal L_{BC|A}=
 \begin{pmatrix}
   k_B-k_A & k_C-k_A\\[2mm]
   h_B-h_A & h_C-h_A
 \end{pmatrix}
 =\begin{pmatrix}
   -1/4&0\\
   \chi+R&2\chi
  \end{pmatrix}.                                  \tag{7}
\]

Therefore

\[
                   \det\mathcal L_{BC|A}=-{\chi\over2}.      \tag{8}
\]

This is independent of `R`.  It is also invariant under adding any common
thermal vector to the three source columns, which is exactly the finite
Schur ambiguity left by `-beta B`.  Thus neither root motion nor slope
motion can restore rank one.

Equivalently, before taking differences the two minors using the
kernel-only column are proportional to `chi+R` and `chi-R`; they cannot
both vanish for nonzero `chi`.  Equation (8) packages that observation in
a root-independent form.

## 4. What this changes, and what it does not

The ordinary four-arm landing space contains at least two independent
channels:

1. **contractible branch merger / canonical-kernel rerouting**, read by
   `D_z g_xy` with no rank jump;
2. **essential-cycle entry versus completion**, read by the opposite signs
   of `D_zE` after the root quotient.

Their exact minor means C4 symmetry and the Schur projection alone do not
force the three-packet leading block to be a pure thermal coordinate.
Consequently Issue #537 should stop the proposed automatic promotion to a
four-packet `R^4 pi_4(R)^4` envelope.  Ordinary four-arm packets already
support the obstruction; an extra-arm estimate cannot be invoked merely
from finite landing algebra.

The certificate is deliberately finite.  It supplies no lower bound on
the probability of these landing types, no assertion that their signed
spatial sum avoids cancellation, no exact-`p_c` to pooled-root transport,
and no continuum field identification.  Any continuation must control the
specific signed ordinary-four-arm functional represented by (7), rather
than replace it by a five/six-arm event.

The exact three-fibre calculation is reproduced by
`scripts/p537_landing_minor_certificate.py`; it imports the already frozen
Bell-eight and torus semantics and performs no enumeration beyond the three
listed configurations.
