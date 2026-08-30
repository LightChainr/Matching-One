# Translation regularity collapses every TM Hall cut

Let `D` be the replicated demand graph and `S` the combined synergy-plus-independent supply multigraph. Torus translations preserve the fixed projective line, layer and transition type, so both weighted graphs are regular.

For every site set `U`, weighted degree counting gives

`g(U)=S_inc(U)-D_in(U)=|U|/N (S_tot-D_tot) + (D_cut(U)+S_cut(U))/2`.

Every term on the right is nonnegative once the all-site aggregate TM cut passes. Therefore the entire induced Hall family follows from one canonical cut. Also

`D_in(U)/S_inc(U) <= D_tot/S_tot`,

so all-site is the worst ratio cut. Positive independent-reservoir edges cross every proper nonempty cut whenever demand is nonzero, making the maximizer unique.

## Uncrossing

The same two-site expansion gives

`g(A)+g(B)-g(A union B)-g(A intersection B)=(D+S)(A\B,B\A)>=0`.

Thus `g` is submodular and the deficit is supermodular. A bisubmodular fallback is unnecessary: translation regularity supplies a stronger exact decomposition.

## Exact bounded census

All 984 rows have regular demand and supply degrees. The oracle checks 2470440 nonempty site cuts, with zero decomposition, Hall or ratio-bound failures. Among 688 nonzero rows, all-site is the unique ratio maximizer every time.
All ratio-`>=9/10` cuts are all-site. The global maximum is `24/25` in 6 N=12 rows. The strongest proper cut is only `576/715` and is a complement-of-one-site cut.

## Revised frontier

The proper-subset Hall family is now a theorem, not an additional conjecture. The remaining topology problem is exactly aggregate TM on arbitrary digital Alexander quotients; connected-set and line-coset reductions add no further burden.
