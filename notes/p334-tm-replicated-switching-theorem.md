# A replicated one-mark switching theorem for TM

Fix a rank-one ambient homology line `ell`, a lower layer `k`, and write
`m=N-k`, `A=|F_k(ell)|`. For a site `v`, let `X_v` count lower states at
which inserting `v` exits from `ell` to rank two.

For an unordered pair `p={v,w}`, `v!=w`, define:

- `D_p`: lower states where both `v` and `w` are exit insertions;
- `Y_p`: lower states where each single insertion stays on `ell`, while the
  double insertion has rank two.

Then `sum_p 2D_p=E_2` and `sum_p 2Y_p=N_new`. The frozen TM inequality is
therefore exactly the cardinality comparison between the following integer
multisets:

- demand: `d_p=2mA D_p` copies of every same-state exit pair `p`;
- synergy supply: `s^Y_p=2mA Y_p` copies of every synergy square `p`;
- independent supply: `(m-1)X_v^2` copies of `(v,v)`, and
  `2(m-1)X_vX_w` copies of `{v,w}` for `v<w`.

Their total difference is

`m A N_new + (m-1)X^2 - m A E_2`.

## The first viable local switch

Requiring a supply token to have exactly the same marked pair as its demand
token is false in 658 of the 984 bounded rows. The first failure is the
matching carrier of `P=diag(2,3)`, `ell=(1,0)`, `k=2`: each of the three
opposite-site pairs has demand 144 but same-pair supply only 24.

The minimal successful relaxation is:

> A demand pair and a supply pair may be matched when they share at least one
> marked site.

Thus a switch preserves one marked insertion and changes only the other.
The fixed `ell` supplies the ambient homology channel; the rule does not use
metric distance, a Gaussian phase, or a quotient-specific fitted ordering.

With lexicographically ordered bins and replica labels, deterministic
shortest-augmenting-path flow gives an explicit integer injection for all 984
rows. It routes 262,842,583 of 466,958,184 demand replicas through a genuine
one-endpoint switch. The most nonlocal-in-pair-space row switches `88/91` of
its demand, showing that the relaxation is structural rather than cosmetic.

## Exact capacitated Hall reduction

Construct a bipartite graph whose left vertices are positive demand pairs and
whose right vertices are typed synergy/reservoir supply pairs. Join two bins
iff their site pairs intersect.

For a family `F` of demand pairs, its supply neighborhood depends only on
`V(F)`, the union of its marked sites. For fixed `V(F)`, adding every positive
demand pair contained in `V(F)` increases left capacity without changing the
neighborhood. Hence capacitated Hall reduces exactly to

`sum_{p subset U} d_p <= sum_{q: q intersects V(D(U))} s_q`

for every site set `U`, where `D(U)` is the positive demand pairs contained in
`U`. This family is necessary and sufficient; integral max flow then lifts to
a token injection by ordering replicas within each bin.

The all-site cut is precisely aggregate TM. Proper subsets are the additional
local content required by a one-mark-preserving switch. Across the existing
atlas, 1,176,258 distinct nonempty induced demand families pass, with no
equality. The tightest ratio is `24/25`, attained by an all-site cut already
known as the tightest aggregate TM row.

## What is now proved, and what remains

Exact theorem: the displayed Hall family is equivalent to existence of the
one-common-site integer injection.

Bounded exact result: every existing HNF line/carrier/layer row satisfies the
Hall family and has an explicit deterministic injection.

Open topological step: prove that digital Alexander rank-one sectors force
all proper-subset Hall cuts for arbitrary quotients. This is strictly sharper
than proving aggregate TM and is now the only remaining switching gap.
