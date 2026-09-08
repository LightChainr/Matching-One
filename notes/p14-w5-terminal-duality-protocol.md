# P14 W5 two-orbit terminal-duality protocol

This protocol fixes the embedding, labels, edge bijections and rejection rules
before producing the committed score certificate.

## Primal disk cell

- Boundary terminals occur counterclockwise as `q = (0,2,1,3)`.
- The internal hub is `h = 4`.
- Rim edge `R_i` joins `q_i` to `q_(i+1 mod 4)` and has probability `r`.
- Spoke edge `S_i` joins `h` to `q_i` and has probability `s`.
- A terminal state is the set partition of the four named primal terminals.
- The exact probability vector is stored in the bivariate Bernstein basis:
  each coefficient counts configurations with `a` open rim edges and `b`
  open spoke edges.

## Two duals that must not be conflated

The four triangular bounded faces are `F_i=(h,q_i,q_(i+1))`; the unbounded
face is `O`.

1. **Spherical dual.**  It has vertices `O,F_0,...,F_3`, with
   `R_i*=(O,F_i)` and `S_i*=(F_(i-1),F_i)`.  Abstractly this is another W5
   and exchanges the two edge orbits.  With the reflection
   `O -> h, F_i -> q_(-i)`, complement-duality is an involution and sends
   `(r,s)` to `(1-s,1-r)`.  Its fixed line is `r+s=1`.

2. **Disk-relative boundary dual.**  The outer face is split into four
   boundary dual terminals `B_i`, one in each arc `(q_i,q_(i+1))`.  It has
   internal vertices `F_i`, boundary edges `R_i*=(F_i,B_i)` and internal
   cycle edges `S_i*=(F_(i-1),F_i)`.  This is the dual whose terminal
   partition is the planar/Kreweras complement of the primal terminal
   partition.

## Frozen tests

The full terminal-vector self-duality gate passes only if the disk-relative
dual is boundary-terminal-preserving isomorphic to the primal W5 and the
configuration-level terminal partition map closes.  A spherical graph
isomorphism alone is not enough.

The certificate will:

- enumerate all `2^8` configurations exactly;
- report the 15 bivariate Bernstein count tables;
- verify disk-relative complementarity configuration by configuration;
- compare degree sequences and boundary-terminal degrees of the primal and
  disk-relative dual;
- verify the spherical complement-duality involution and weight preservation
  on `r+s=1`, then test whether its output partition is a function of the
  input terminal partition;
- if it is not, choose the smallest obstruction by minimizing the maximum
  open-edge count of a witness pair and then lexicographically sorting bit
  strings in the fixed edge order `R_0,...,R_3,S_0,...,S_3`.

No equality such as `P(0123)=P(0|1|2|3)` is assumed.  No scalar root from the
parent screen is interpreted as a threshold, periodic comparison or bound.
