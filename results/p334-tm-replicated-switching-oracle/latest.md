# Replicated one-mark switching for TM

Write the TM determinant as a doubled-unordered supply-demand problem. A same-state exit pair `{v,w}` has `2mA` demand replicas. A synergy square has `2mA` supply replicas. The independent reservoir contains `(m-1)` replicas of every ordered pair of exit marks. Its total capacity is `(m-1)X^2`.

## The local rule found by canonical flow

Locking supply to exactly the same site pair fails. The minimal bounded failure is `N=6`, matrix `[[2, 0], [0, 3]]`, matching carrier, line `[1, 0]`, lower layer `2`. Overall, exact-pair locking fails on 658 of 984 rows.

The first successful locality is one-common-site compatibility: preserve either `v` or `w` and switch only the other marked site. Lex-first integral max flow succeeds on every row. This rule uses the fixed ambient line `ell`, but no metric distance, quotient-specific phase, or fitted parameter.

## Exact Hall family

For a set `U` of marked sites, let `D(U)` be demand on every positive demand pair contained in `U`. Let `S(U)` be supply on every synergy or reservoir pair incident to a site actually used by `D(U)`. Then

`D(U) <= S(U) for every U`

is necessary and sufficient for the one-common-site injection. Necessity is Hall. For sufficiency, the neighborhood of any demand family depends only on its union of marked sites, and adding every positive demand pair inside that union maximizes left capacity without changing its neighborhood. Integrality of bipartite max flow gives a token injection.

The oracle checks 1176258 distinct nonempty induced demand families with zero failures. The tightest cut has `required/available=24/25` and is the all-site TM cut. The largest canonical one-endpoint switching fraction is `88/91`.

## Boundary of the theorem

This closes the switching problem on the bounded atlas and replaces an opaque aggregate inequality by an exact local Hall criterion. What is not yet proved is that digital Alexander topology forces every proper-subset Hall cut on arbitrary quotients. That is now the sole topological gap; no larger-N evidence is needed to state it precisely.
