# Leaving the closed-source endpoint through one saturation defect

**New finite mechanism.** Deleting a saturated site changes the closed action
by `Delta S=3−2 k_null−ell`: k_null is the number of lost zero-ambient-image
graph cycles and ell the ambient-rank loss. At this endpoint ell is at most1;
only an alternating child face can produce it. The induced child defect is
a change of local diagonal connection, not an ordinary independently occupied
child site.

We retain `S=C+F+Bv`, the exact action and nonalias square-torus convention in
[the closed-source proof](checkerboard-positive-source-closure.md). The
[b8d043fc overview](https://github.com/LightChainr/Matching-One/blob/b8d043fc/notes/decimation-closed-source-and-global-u.md)
already supplies the normalized one-hole insertion formula. Here its weight,
cycle/rank meaning and fixed-t thermal derivative are made explicit. No new
numerical outcome, enumeration or sampling is used.

## 1. Exact loss of cycles and winding

Let X_plus have all M=N/2 checkerboard A sites occupied; X_minus_a deletes
one A site a. Before deletion let d be its occupied NN degree. Let c be the
number of components into which its old occupied component splits after
removing a; set c=0 for an isolated a. Every surviving component contains
one of its former occupied neighbors, so

```text
Delta K=−1,    Delta Bocc=−d,    Delta C_B=c−1,
beta1=Bocc−K+C_B,    k:=beta1_plus−beta1_minus=d−c.
```

The remaining graph's cycle space is a subspace of the original one. Set
`ell=r_plus−r_minus=−Delta q`. Then `0≤ell≤k` and the loss of the kernel of
the ambient-H1 map is exactly `k_null=k−ell≥0`. Using the already proved
`S=2 beta1−3K−q+2N` gives

```text
Delta S = 5−2d+2 Delta C_B−Delta q
        = 3−2k+ell = 3−2k_null−ell,
q_minus=q_plus−ell,
E_minus−E_plus=ell²−2ell*q_plus,
exp(t Delta S)=exp(3t) exp(−2t k_null) exp(−t ell).
```

Thus a lost null cycle and a lost ambient cycle have different exact source
costs. The constant exp(3t) does not remove the other factors or the t-dependent
endpoint distribution. This is not a proof of one overall response gain.

## 2. Which single holes can change rank at saturation?

The four B-neighbors of a are the four corners of one child face. If two
occupied B-neighbors are adjacent around a, their two-edge path through a
has an alternative path through the other, still-filled A corner of the
original unit square. These paths have the same lifted displacement.
For three or four occupied B-neighbors the corresponding boundary paths
connect all of them, so removing a changes only zero-ambient cycles.
Zero or one occupied neighbor cannot support an ambient cycle through a.

The sole exceptional pattern has exactly two **opposite** occupied
B-neighbors. In the complemented child its four corners alternate occupied
and vacant. Removing a deletes the matching diagonal between the two
vacant child corners and inserts the complementary diagonal between its
two occupied corners. This is a local change of connection convention.
It leaves all child site labels unchanged. The occupied parent loses at
most one cycle, so ell is0 or1, never2.

The complete possibilities are summarized without enumerating configurations:

| Occupied neighbors of a | c | k | ell | Delta S |
|---|---:|---:|---:|---:|
| d=0, isolated |0|0|0|3|
| d=1 |1|0|0|3|
| d=2, adjacent |1|1|0|1|
| d=2, opposite, separates the occupied component |2|0|0|3|
| d=2, opposite, component remains connected |1|1|0 or1|1+ell|
| d=3 |1|2|0|−1|
| d=4 |1|3|0|−3|

In the connected opposite-pair case, ell distinguishes an ambient-redundant
cycle from loss of an independent winding direction. Determining that case
requires the exterior connection and its gain: degree or the four occupancy
bits alone do not determine it. The new diagonal is only a topological
description; the exact source weight also includes the cycle and component
changes above, including rank-preserving holes.

Rank loss can occur without a degenerate quotient. On the child L×L axis
torus, L≥5, occupy exactly
`U={(0,0),(2,0),(3,0),...,(L−1,0),(2,1),(1,1)}`.
This is one NN path of ambient rank0. Adding the diagonal(0,0)–(1,1)
in its alternating face closes a cycle of winding(−1,0), giving rank1.
For parent periods `(1+i)L Z[i]`, deleting A site a=1+i performs exactly
that switch. The parent rank drops2→1, with k=ell=1 and Delta S=2.
This is a paper configuration, not a simulated example.

## 3. Normalized finite-t insertion, including rank-preserving weights

For0<p<1 use `epsilon=1−s`, `p_A=1−epsilon(1−p)`, `p_B=p` and fixed real t. Let
`<.>_+` denote the normalized saturated law with weight exp(t S_plus),
`mu_O=<O_plus>_+`, and `w_a=exp(t Delta_a S)`. The partition and observable
numerator expansions relative to the saturated partition function are

```text
Z_epsilon/Z_+ = 1+epsilon(1−p)[sum_a <w_a>_+−M]+O(epsilon²),
Z_epsilon[O]/Z_+
 = mu_O+epsilon(1−p)[sum_a <w_a O_minus_a>_+−M mu_O]+O(epsilon²).
```

Their ratio recovers the overview's exact derivative

```text
j_O:=partial_epsilon <O>|0 = (1−p) h_O,
h_O = sum_a <w_a(O_minus_a−mu_O)>_+
    = sum_a [<w_a Delta_a O>_+ + Cov_+(w_a,O_plus)].
```

The second term must be retained. Even when ell=0 and q/E do not change
on that deletion, its source weight can be correlated with the global
observer. At t=0 all w_a=1, this covariance term vanishes, and only genuine
rank-changing alternating faces contribute to q/E. In particular

```text
j_q(p,0)=−(1−p) sum_a <ell_a>_+ ≤0,
j_E(p,0)=(1−p) sum_a <ell_a²−2ell_a q_plus>_+.
```

This sign for the fixed-p matching-mean response at t=0
does not determine the sign of the root-comoving U response or of the
fully normalized response at arbitrary t.

## 4. Fixed-t thermal derivative and the original U

Only B-site occupations vary with p at the endpoint. Put K_B for their
count and u=p(1−p). Since w_a and the configuration observables have no
explicit p dependence at fixed t,

```text
h'_O = (1/u) sum_a [Cov_+(w_a(O_minus_a−mu_O),K_B)
                    −<w_a>_+ Cov_+(O_plus,K_B)],
j'_O = −h_O+(1−p)h'_O.
```

The minus h_O term is forced by the defect dose1−p. It also accounts for
the p_A derivative induced when differentiating the declared interior path;
it cannot be removed by treating the hole probability as a constant dose.

Use the original fixed parent geometry weights, `Q=mean(q)`, `Y=P4(E)`,
and on a simple saturated pooled-root branch put `D=Q_p`, `r=Y_p/D`,
`A_N=N^(13/8)/2`. Construct j_Q,j_Y and h_Q,h_Y linearly from the formulas
above, with the same t and endpoint population. Then

```text
partial_epsilon p0 = −j_Q/D,
V_epsilon(t):=partial_epsilon U|0
 = A_N/D * partial_p(j_Y−r j_Q)
 = A_N/D * [(1−p)partial_p(h_Y−r h_Q)−(h_Y−r h_Q)] at p0(t).
```

The derivative acts on r as well. The endpoint s derivative has the opposite
sign, `partial_s U|1=−V_epsilon`. These formulas include root relocation,
thermal-slope normalization and the p-dependent defect dose.

The defect can be represented exactly on a modified child connection
structure with the displayed weights, but it is not the unchanged ordinary
child measure with a shifted independent site probability. Endpoint closure
alone imposes no thermal-only relation for h_Q,h_Y and no value of its mixed
t derivative. Whether the resulting global-U response reduces to a common
thermal tangent or a source-independent gain remains a distinct test; no
numerical value, field identification or additional descriptor is supplied
here.
