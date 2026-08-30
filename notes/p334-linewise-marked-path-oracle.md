# Ambient-H1 marked paths behind BA and TM

Fix a projective primitive line `ell` in ambient `H1`. A birth boundary
`0->ell` selects this line. An exit boundary `ell->H1` creates the unique
nonzero quotient direction in `H1/ell`; no arbitrary transverse basis is
needed.

## BA: joint paths versus independent marks

On one fixed-line layer, let `b(T)` and `x(T)` be the birth and exit degrees.
Then

`C=sum_T b(T)x(T)`

counts jointly marked two-step forks `0<-T->2`. BA is

`A C >= B X`, where `B=sum b`, `X=sum x`.

Equivalently,

`2(AC-BX)=sum_(S,T)(b(S)-b(T))(x(S)-x(T))`.

The positive terms are concordant marked state pairs and the negative terms
are discordant pairs. Thus BA asks for a weight-preserving injection from
discordance to concordance. It cannot be replaced by pointwise comonotonicity,
which already fails at `N=8`.

Alexander complement gives one exact switch. If `T` carries primal marks
`(birth v, exit w)`, then `E\T` carries matching marks `(birth w, exit v)` at
the reflected layer. This preserves `C`, swaps `B,X`, and preserves the BA
determinant. It shows the two carrier BA problems are the same marked problem,
but does not yet construct the within-layer discordance injection.

The bounded atlas contains eight nonzero BA equalities. In each,
concordance and discordance masses are both 1600. Any general proof must allow
a genuine bijection in those cases; a strict-surplus argument cannot work.

## TM: new transverse squares and independent exit pairs

Set `m=N-k`. For a lower state `S`, an internal insertion `v` may make a second
site `w` newly exit-pivotal. Then

- `S+v` stays on `ell`;
- `S+w` also stays on `ell`;
- `S+v+w` has rank two.

This is a synergy square, and `(S,v,w)<->(S,w,v)` is a canonical involution.
Let `N_new` count its oriented sides. Let

`E_2=sum_S x(S)(x(S)-1)`

count ordered distinct exit marks at one lower state. Direct expansion of the
TM determinant gives the exact path inequality

`m A N_new + (m-1)X^2 >= m A E_2`.

The first supply is replicated synergy squares. The second is a replicated
reservoir of independent exit-mark pairs. Both are essential: 74 audited rows
have `N_new=0` but `E_2>0`. The six tightest `N=12` rows also have no synergy
contribution; independent supply 43200 covers demand 41472, a ratio `24/25`.
Therefore a synergy-only injection is false even though the square involution
itself is exact.

## Weakest verified topology theorem

For every primitive line `ell`, require:

1. concordance mass is at least discordance mass;
2. replicated synergy-plus-independent TM supply is at least ordered-exit
   demand;
3. the same conditions hold on the complement carrier.

These are local, integer and directly auditable. By the two-carrier moment
theorem they imply ULC. All 984 existing line/carrier/layer rows satisfy them,
all 59,922 complement state marks swap exactly, and the synergy-square
involution has zero failures.

What remains open is a topology-derived switching for the replicated marked
sets. The exact obstructions show its necessary shape: it must be aggregate,
permit BA bijections, and use the independent `X^2` reservoir in TM.
