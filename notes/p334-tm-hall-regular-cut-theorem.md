# The regular two-mark cut theorem

The replicated TM injection of `551f890` produced a capacitated Hall family
indexed by site sets. That family collapses to the ordinary all-site TM
inequality because finite-torus translation symmetry makes both pair graphs
regular.

## Setup

Let `V` contain `N` sites. Let `D` be the nonnegative weighted graph of
replicated same-state exit-pair demand, and let `S` be the nonnegative weighted
multigraph of synergy plus independent-reservoir supply. Supply loops are
allowed and count twice in weighted degree.

For `U subset V`, write:

- `D_in(U)` for demand pairs with both ends in `U`;
- `S_inc(U)` for supply pairs with at least one end in `U`;
- `D_cut(U),S_cut(U)` for pair weight crossing `U,V\U`.

The finite quotient translation group acts transitively on sites. Translation
preserves the carrier, subset size, fixed projective line `ell`, and every
rank-transition type used to define `D` and `S`. Consequently both weighted
graphs are regular.

## Exact cut decomposition

Weighted degree counting gives

`2D_in(U)+D_cut(U)=2|U|D_tot/N`

and

`S_inc(U)=|U|S_tot/N+S_cut(U)/2`.

Subtracting proves the exact identity

`g(U)=S_inc(U)-D_in(U)`

`    = |U|/N (S_tot-D_tot) + (D_cut(U)+S_cut(U))/2`.

Therefore `D_tot<=S_tot`, which is precisely aggregate TM, implies every Hall
cut. It also gives

`D_in(U)/S_inc(U) <= D_tot/S_tot`

whenever demand is nonzero. If demand is nonzero then translation regularity
makes every `X_v` positive, and `m>=2`; hence the independent reservoir has a
positive off-diagonal edge across every proper nonempty cut. The ratio
inequality is then strict for proper `U`. Thus all-site is the unique worst
ratio cut.

This is stronger than reductions to connected or line-coset sets: those
classifications are unnecessary.

## Uncrossing identity

For arbitrary `A,B subset V`, a two-endpoint truth table gives

`g(A)+g(B)-g(A union B)-g(A intersection B)`

`    = (D+S)(A\B,B\A) >= 0`.

Thus `g` is submodular and the Hall deficit is supermodular. Dangerous sets
uncross under union/intersection, but regularity goes further by identifying
the unique canonical maximizer directly.

## Bounded exact classification

Across all 984 existing rows:

- both weighted degree vectors are regular in every row;
- 2,470,440 nonempty site cuts satisfy the exact decomposition;
- 688 rows have positive demand and all-site is the unique ratio maximizer in
  every one; the remaining 296 rows have zero demand;
- every cut with ratio at least `9/10` is all-site;
- the global maximum `24/25` occurs in six `N=12` primal rows;
- the strongest proper cut is `576/715`, obtained by deleting one site from a
  tight `N=12` geometry.

## Revised theorem frontier

The induced Hall family is no longer an independent conjecture. For every
translation-homogeneous torus row it is an exact corollary of aggregate TM.
The only remaining TM problem is to prove the all-site two-carrier moment
inequality for arbitrary digital Alexander quotients.
