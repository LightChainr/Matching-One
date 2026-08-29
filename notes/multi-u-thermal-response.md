# Multi-u thermal-response templates

Status: C0 definition freeze for issue #119, C1 N=10 oracle. Descriptive P49 only. Not a P43 / Issue #57 target.

## Why freeze the shape before the next full curve

Centre derivatives confound ordinary analytic, Jordan/log, and nonlinear coordinate corrections. The same simulations already solve intrinsic levels `u ∈ {0, 0.025, 0.05}`. Issue #101 isolates the even nonlinear thermal coordinate as

```text
Q_N = c_{0.05} - c_{0.025} ~ N^{-3/4}
```

and the odd width as `w_u ~ u N^{-3/8}`. Across the frozen grid this is the two-component template

```text
w_u / u           ~ B N^{-3/8}
(c_u - c_0) / u^2 ~ A N^{-3/4}
```

Ordinary q=2 analytic corrections add higher odd/even powers of `u`. A rank-2 Jordan factor multiplies the same u-shape by `log N` and therefore changes the N-power of the even piece.

Do not add levels after looking. The vector over `u` is one correlated observation: recompute `p_±^u` inside each delete-one replicate.

## Exact oracle

On the C4 N=10 Beta(3,3) control every midpoint is `1/2`, so the even u-shape vanishes identically. That is the odd-around-1/2 check of the template, not a continuum amplitude.

## Descriptive P49

N=130 and N=170 are children of different doubling lineages. They are used only to see that `w_u/u` and `(c_u-c_0)/u^2` are already u-stable at current sizes. They are not a doubling test and not a P43 score.

A later 145→290 full curve may use the frozen monomials as a joint template, provided the level solver is recomputed inside each jackknife replicate.
