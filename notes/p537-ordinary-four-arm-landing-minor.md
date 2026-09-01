# Ordinary four-arm landings have two raw jump symbols; the Schur/P4 minor remains open

## Corrected decision

Three exact `4 x 4` NN-torus fibres show that ordinary four-arm connectivity
supports two independent **raw jump symbols**:

1. canonical-kernel rerouting with no ambient-rank jump;
2. rank entry/completion with no kernel jump.

If `chi` is inserted as a formal coefficient of `E`, their raw difference
matrix has determinant

\[
                              -\chi/2.                         \tag{1}
\]

This proves that arm number and the local C4 action alone do not make the
kernel and readout channels proportional.

It does **not** yet give the Schur- and P4-projected minor requested by Issue
#537.  The earlier version conflated two interfaces:

* equal background `K` makes `S_i,B_i` identical on the three two-state
  fibres, but `-beta B_i` is contracted with the fibre-dependent `H_i`, so
  it is not a common additive source coordinate;
* a C4 orbit inside one square quotient is not the axis-minus-tilted P4
  projector.  A nonzero `chi` must come from paired, separately normalized
  geometry amplitudes.

The shortest proved statement is therefore:

> A rank-one result, if true, must be a nontrivial identity of the complete
> Schur bilinear and paired geometry transfer.  It cannot follow from
> ordinary-four-arm counting plus local C4 symmetry.

No claim is made that an extra fifth/sixth arm is unnecessary after the full
projection.

## 1. Raw landing coordinates

For a fixed pair `(x,y)` and thermal site `z`, define

\[
 k=D_zg_{xy},\qquad h_R=D_z(y-Rm)=\chi D_zE-RD_zq.            \tag{2}
\]

The first coordinate is the kernel-reconnection jump.  The second is the
formal readout jump after the row operation `y -> y-Rm`.  Entry and
completion obey

\[
\begin{array}{c|cc|c}
 &D_zq&D_zE&h_R\\ \hline
 0\to1&1&-1&-\chi-R\\
 1\to2&1&+1&+\chi-R.
\end{array}                                                   \tag{3}
\]

A rank-preserving reconnection has `D_zq=D_zE=0` but may have `k!=0`.
These are the two exact midpoint-derivative channels, but the pair of jumps
is not itself the full mixed Hessian.

## 2. Three physical ordinary-four-arm fibres

All coordinates are modulo four.  Unlisted sites and `x,y,z` are vacant.
Each background has `N=16`, `K=6`.  At every `z`, the incident
`(N,E,S,W)` occupation is `(1,0,1,0)` or `(0,1,0,1)`: exactly four
alternating incident arms and no local fifth branch.

### Entry `A`

```text
occupied = {(0,0),(1,0),(0,1),(2,0),(1,1),(0,2)}
z=(3,0),       (x,y)=((2,1),(1,2))
rank: 0 -> 1,  16 g_xy: 0 -> 0.
```

Hence `(k,h_R)_A=(0,-chi-R)`.

### Kernel rerouting `B`

```text
occupied = {(0,0),(1,0),(0,1),(2,0),(1,1),(2,2)}
z=(2,3),       (x,y)=((2,1),(1,2))
rank: 0 -> 0,  16 g_xy: 4 -> 0.
```

Hence `(k,h_R)_B=(-1/4,0)`.

### Completion `C`

```text
occupied = {(0,0),(1,0),(0,1),(2,0),(0,2),(3,0)}
z=(0,3),       (x,y)=((1,1),(2,1))
rank: 1 -> 2,  16 g_xy: 4 -> 4.
```

Hence `(k,h_R)_C=(0,chi-R)`.

The distinction is in the landing identification of the two occupied arms:
closing the first essential cycle, merging two zero-image branches, or
closing the second essential cycle.  An arm count alone does not identify
these outcomes.

## 3. What local C4 proves

Rotate each complete tuple `(occupied,x,y,z)` through its C4 orbit.  Ambient
rank, canonical `Kreg` and entry/completion type are unchanged, so the raw
values above survive local orbit averaging.

Subtract fibre `A` from `B` and `C` at the level of the jumps (2).  This gives

\[
 \mathcal L^{raw}_{BC|A}
 =\begin{pmatrix}
   -1/4&0\\
   \chi+R&2\chi
  \end{pmatrix},
 \qquad \det\mathcal L^{raw}_{BC|A}=-{\chi\over2}.            \tag{4}
\]

This is independent of `R`, but it is only a formal raw-symbol minor.  A C4
quarter-turn has spin-four character `+1`; that observation does not supply
the nonzero axis-minus-tilted coefficient `chi`.

## 4. Why common `K` does not remove the Schur term

The exact projected quantity is

\[
 T_t=\left\langle H,(a-Ea)S-\beta B\right\rangle_{pool}.      \tag{5}
\]

For a two-state fibre with off-`z` occupation `K=6`, the values
`S_0,S_1,B_0,B_1` are indeed common across `A,B,C`.  But its conditional
Schur contribution is

\[
 (1-p)H_0[-\beta B_0]+pH_1[-\beta B_1],                      \tag{6}
\]

and `(H_0,H_1)` differs among entry, rerouting and completion.  Therefore
column subtraction does not cancel (6).  Common `K` controls the thermal
score, not the bilinear contraction.

The same problem affects centering: `Ea` is global and geometry-specific,
while `a_0,a_1` differ among the fibres.  The jump `D_za` alone does not
retain these midpoint terms.

## 5. The actual finite matrix still required

For each landing state `ell` and geometry `g`, retain the complete tensor

\[
 (q_0,q_1,E_0,E_1,a_0,a_1,K_{-z})_{g,\ell},                 \tag{7}
\]

its outer landing amplitude, and the geometry-specific means entering
`H`, `Ea`, `R`, and `beta`.  Form

\[
 \sum_{i=0}^1p_iH_i\{(a_i-Ea)S_i-\beta B_i\}                \tag{8}
\]

before separately normalizing axis and tilted geometries and applying P4.
A transfer minor needs two source boundary states and two thermal/readout
boundary states; three already-contracted whole fibres are not automatically
that matrix.

The current certificate remains useful as a minimal obstruction candidate:
ordinary four-arm landings contain both raw channels.  It supplies no
projected nonzero minor, spatial probability bound, signed-sum result,
exact-`p_c` transport, or continuum identification.

The exact raw calculation is reproduced by
`scripts/p537_landing_minor_certificate.py`; it checks only these three
configurations.
