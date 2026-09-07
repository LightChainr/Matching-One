# Experiment-complexity atlas — probe round 1 slice

A deliberately incomplete first slice of the atlas requested by #599 (direction
1).  Every number below is tied to an exact finite computation in the
repository or to this probe's verified scripts.  "State dimension" is not a
column: it is what this atlas is testing.

Evidence labels: **P** = proved exact (repository), **p** = proved exact (this
probe, files cited), **n** = numerical/empirical at fixed width, **O** = open,
**C** = claim/counterexample boundary, **F** = false as stated.

## 1. P398 (widths 4–8 unless noted)

| notion | w4 | w5 | w6 | w7 | w8 | anchor | evidence |
|---|---|---|---|---|---|---|---|
| microscopic states `n_micro` | 14 | 42 | 132 | 429 | 1430 | Catalan(w) | P |
| C2-orbit quotient `n_sym` | 10 | 26 | 76 | 232 | 750 | `(Cat(w)+C(w,floor w/2))/2` | P (#597/#598, re-verified Gate-1) |
| exact positive / strong-lumping `r_pos` | 10 | 26 | 76 | 232 | 750 | coarsest admissible lumping **=** orbit partition block-for-block | P (Gate-1 branch A) |
| signed linear realization `r_lin` | 10 | 26 | 72 | 218 | ≥150 | exact mod-p elimination | P/O (w8 open, #593) |
| Krylov transport `r_transport` | 4 | 6 | 8 | 8 | 8 | #580/#594 | n |
| balanced order (frozen tol) | 3 | 3 | 4 | 4 | 4 | #594; **full = quotient** | n + p |
| memory kernel effective order | 2 | 3 | 4 | 4 | 4 | 99.9% energy (numerical order 4,9,12,13,14) | n (#588) |
| task singular spectra (this probe) | full ≡ quotient to 1e-16 | | | | | probe note | p |
| `halves_linked` parity | even | **odd** | even | **odd** | even | C2 label | P |

Recorded non-implications (P398):

* `r_pos < r_lin`?  False at w4–7 (`r_lin <= r_pos`), direction flips open at
  w8.  So even *within* linear-algebraic notions the order is not fixed.
* `r_lin <= r_pos` holds w4–7; `r_lin` is not an orbit count (72 ≠ 76 at w6).
* balanced order is **not** a function of `n_micro` nor of `r_pos`; it is a
  function of the (source, readout, horizon) task after the exact quotient
  factor (probe note; equality full ≡ quotient to 1e-16).
* symmetry quotient is an **exact first factor**; balancing is a second,
  task-relative compression.  The two claims are not in competition.

## 2. Cut-network / predictive language

| notion | value | anchor | evidence |
|---|---|---|---|
| complete unbranched survival law `S(z)` | identical across #435/#549 bank, independent of `a` | `S_k(z)=S(z)^k` | P |
| single-round fork class count | `k+1` classes (`a=0..k`) inside one survival class | `F_{k,a}` closed form | P |
| fork response rank (any finite single-round language) | `<= 2` while classes `= k+1` | probe note (`affine in a`) | p |
| response rank with polynomial degree-d tests | `<= d+1` | Vandermonde lemma | p |
| r=1 bounded summaries vs depth-2 compositional language | insufficient (witness gap 1/525) | #550 | P |
| minimal depth for unbounded rank | **O** | equals depth producing degree-(k+1) probability | O |
| `q3` in `x` (one-step coordinate) | degree 3 → rank capacity 4 | #401 formula | P |
| Nerode-class count vs Hankel rank vs nonnegative rank | separate; classes can exceed rank | probe note | p |

## 3. Finite-size flow and threshold claims

| claim | status | anchor |
|---|---|---|
| 99% one-direction quantile flow ⇒ 1-d RG state | **F** (generic smooth tangent) | #582/#584 null, accepted |
| bounded finite-horizon task order ⇒ threshold constraint | **F** in general; bridge needed | probe note (unbranched ≡ but fork ≠) |
| `E_p[X]=0` homological balance implied by low task order | **O** / no known implication | probe note |
| structured remainder reproducible by smooth family | **O** — proposed discriminator (cocycle/curvature/crossed labels) | probe note |

## 4. Known unknowns ranked by cheapest falsification

1. **w8 exact `r_lin`** (is it > 218, does `r_lin > r_pos` happen?) — one exact
   mod-p elimination at n=1430, cheap.
2. **two-round #549 fork probability degree in `a`** — one enumeration in the
   protocol, decides rank-3 reachability.
3. **halves_linked as the unique parity-odd dictionary element at other
   dictionaries D1/D2** — direct check, clarifies "marked readout" reachability.
4. **w9/w10 memory order saturation vs slow growth** — needs the 4862/16796
   state computation (documented heavy; #593).
