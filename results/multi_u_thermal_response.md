# Multi-u thermal-response templates

Source: `scripts/multi_u_thermal_response.py`.
Claim level: C0 template freeze, C1 N=10 oracle. P49 rows are development-only.
Not a P43 / Issue #57 target. Frozen `u={0, 0.025, 0.05}`; do not add levels.

## Frozen monomials

```text
w_u / u           ~ B N^{-3/8}
(c_u - c_0) / u^2 ~ A N^{-3/4}
```

These are the leading #101 coordinate-nonlinearity shapes. Ordinary q=2
analytic corrections add higher powers of u. A Jordan log multiplies the
same u-shape by log N. Templates must be used jointly with delete-one
covariance; they are not three independent tests.

## N=10 Beta(3,3) oracle

Even midpoints vanish: `Q_10=0.0`, every `c_u=1/2`.

## Descriptive P49 N=130/170

These sizes are children of different doubling lineages, not a doubling pair.

| N | `w/u` at 0.025 | `w/u` at 0.05 | `(c-c0)/u^2` at 0.025 | `(c-c0)/u^2` at 0.05 |
|---:|---:|---:|---:|---:|
| 130 | 0.09224507 | 0.09229116 | -0.002149 | -0.002151 |
| 170 | 0.08347380 | 0.08351584 | -0.001787 | -0.001789 |

On both sizes `w_{0.025}/w_{0.05}≈1/2` and the two even ratios agree to
relative O(10^{-3}). That is compatible with the frozen monomials. It is
not a prospective Jordan discriminator and is not a P43 score.
