# Certified finite-size bound feasibility gate

## Decision

**Methodologically interesting but weak.** The Riordan--Walters reduction applies directly
to square-site percolation, and its binomial certificate is inexpensive. The dominant cost is
making the block event sufficiently likely at a parameter close to `pc`; the repository does
not yet implement this open-boundary event, and a width `5e-4` interval already exists.

## Certifiable theorem shape

Partition the plane into `s x s` cells. A renormalized bond is open when both adjacent cells
have unique largest open clusters and those clusters connect in their `2s x s` union.
Nonincident block bonds depend on disjoint site sets, so the block model is 1-independent.
If its edge probability is strictly above `0.8457`, the modern deterministic bound on
`p_max(Z^2)` implies percolation of the block model and hence of the original model.

For `N` independent trials and `M` successes, reject `H0: pi<=0.8457` using the exact tail
`P[Bin(N,0.8457)>=M]`. With family-wise error `1e-6`, three attempts on each of the original
and matching lattices give per-run alpha `1/6000000`.

## Exact statistical thresholds

| protocol | N | required successes | null tail |
|---|---:|---:|---:|
| 2007 constant `0.8639` | 400 | 378 | `1.14899035289400922E-7` |
| 2022 constant `0.8457` | 400 | 373 | `9.51451513592604452E-8` |

The first row exactly reproduces the published `400/378` test. The improved deterministic
constant lowers the cutoff to `373`, but 400 trials have useful power only when the true block
event probability is already high (about the mid-0.9 range).

## Two-sided interval

A successful square-site run at `p` gives `pc<=p`. A successful run on the square matching
site graph at `q` gives `pc(square)>=1-q`. These are high-confidence statistical bounds, not
deterministic inequalities, and require predeclared tests plus independent genuine Bernoulli trials.

## Heuristic cost warning

Using the paper's non-rigorous thermal-window scaling `delta p ~ s^(-3/4)`, reducing interval
width by a factor `r` multiplies linear scale by `r^(4/3)` and area work by `r^(8/3)`.
This estimate is for planning only and is not part of the certificate.

| target width | linear-scale multiplier | area-work multiplier |
|---:|---:|---:|
| `5.0e-04` | `1.000` | `1.000` |
| `2.5e-04` | `2.520` | `6.350` |
| `1.0e-04` | `8.550` | `73.100` |
| `5.0e-05` | `21.544` | `464.159` |
| `1.0e-05` | `184.202` | `33930.220` |

## Production gate

Implement and independently validate the exact open-boundary block event first. Exploratory
samples may choose `s,p`, but final certification samples must be fresh. Stream separation and
reproducibility are necessary software controls; they do not turn a pseudorandom generator into
a mathematical source of genuine randomness.

Sources: Riordan--Walters `arXiv:math/0702232`; Balister--Johnston--Savery--Scott
`arXiv:2206.12335`.
