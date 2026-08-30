# Linewise marked-path form of BA and TM

Fix a primitive ambient homology line `ell`. A birth mark is a boundary `0->ell`; an exit mark is a boundary `ell->H1`, equivalently creation of the unique nonzero quotient direction in `H1/ell`.

## BA as concordance of marked two-step paths

At one layer, `sum b(T)x(T)` counts joint marked paths `0<-T->2`. The exact determinant

`A sum bx - (sum b)(sum x)`

is half the concordance mass minus discordance mass over ordered state pairs. BA is precisely the statement that concordant boundary fragilities dominate discordant ones. Alexander complement switches `(birth v, exit w)` to `(birth w, exit v)` on the matching carrier and preserves this determinant.

## TM as synergy squares plus an independent reservoir

Let `m=N-k`, `X=sum x(S)`, `E_2=sum x(S)(x(S)-1)`, and `N_new` count oriented triples `(S,v,w)` in which both single insertions preserve `ell` but the double insertion creates rank two. Then TM is exactly

`m A N_new + (m-1) X^2 >= m A E_2`.

The canonical switch `(S,v,w)<->(S,w,v)` proves every new-exit event is an oriented side of a synergy square. It does not prove TM by itself: the independent `X^2` reservoir is essential.

## Exact bounded audit

All 984 line/carrier/layer rows pass both path inequalities. The square involution has 0 failures, and complement swaps birth/exit degrees on 59922 states with 0 failures.
BA reaches the nontrivial equality `discordance/concordance=1` in 8 complement-paired rows, so a proof must allow bijection rather than strict surplus.
TM is nearly tight: maximum demand/supply is 24/25 in 6 N=12 rows. In those maximizers the synergy contribution vanishes. More broadly, 74 rows have no new synergy square but positive same-state exit-pair demand.

## The local topology theorem and the remaining switch

If the BA concordance inequality and the replicated TM path inequality hold for each primitive line on both complementary carriers, the two-carrier moment theorem gives ULC. This is the weakest currently verified linewise condition. A global statewise injection is too strong; the open topological work is a weight-preserving switching of discordant BA marks and a replicated TM injection that uses both synergy squares and independent exit pairs.
