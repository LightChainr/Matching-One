# Root-count algebra of the #549 fork language: rank growth needs depth, not width

**C2 note (Program B of #599).**  An exact mechanism analysis of the #549
parallel-gadget fork family, answering #599 Program B Q1/Q2/Q3 at the level of
the published protocol, plus a small exact abstract model that isolates the
mechanism and a precise statement of what would be needed to lift it to the
real N16 network.

All arithmetic is exact (`Fraction`); no sampling.  Script:
`scripts/probe/c2_root_rank.py`.

---

## 1. Protocol reconstruction (from the #549 branch verifier)

The #549 verifier builds the fork probability from three exact constants of
the #435 N16 base gadget and two per-gadget successor tables:

```
BASE_SAFE = (1,7,18,20,8,0,0,0,0)          # safe m-subset counts (one gadget, 8 future vertices)
EXIT_A    = (1,1,1,2,2,3,3)                # safe-successor exit counts, A-type gadget
EXIT_B    = (1,2,2,2,2,2,2)                # safe-successor exit counts, B-type gadget
```

with `sum EXIT_A = sum EXIT_B = 13` (equal first successor moment), while
`sum x^2` is 29 (A) vs 25 (B) (second moment distinguishes the types).  The
protocol count is

```
F_{k,a} = [ sum over gadgets g of sum over vertices v in g of (7k - x_v)^2 ]
          / [ 8k (8k-1)^2 ],
```

which is the successor-sum form of: pick **one shared-update root vertex**
uniformly (8k choices); then two clones each independently pick a safe
second vertex (each `7k - x` safe choices out of `8k-1`); success if both are
safe.  This reproduces the closed form

```
F_{k,a} = [343 k^3 - 182 k^2 + 25 k + 4a] / [8k (8k-1)^2],
gap = F_{k,a+1} - F_{k,a} = 1/[2k(8k-1)^2],
```

for every `k = 1..12`, `a = 0..k` (re-verified here exactly).

**Reading:** the protocol has exactly **one root**, whose type (the gadget it
lands in) is the only random variable the whole event conditions on.  The two
clones are conditional on that single root.

## 2. Why the single-root fork language caps rank at 2

**Proposition (single-root ⇒ affine in the class coordinate).**  In a parallel
bank of `k` gadgets whose class is the count `a` of A-type gadgets, any test of
the form "one shared-update root, then any number `q` of clones whose success
is conditional on that root and factorizes over gadgets" has success
probability affine in `a`.

*Reason.*  Conditional on the root falling in a gadget of type `tau`, the
clone-success factors separate over gadgets and each gadget contributes a
factor that depends only on its own type.  Summing over the root position
gives `a * C_A + (k-a) * C_B` plus a term independent of the type layout, i.e.
an affine function of `a`.  The clone count `q` does not change this: more
clones still condition on the same single root (`#549` itself has `q = 2` and
is affine; the successor-sum form shows the `q`-clone variant would sum
`s_v^q` per root vertex, still gadget-wise additive, hence still affine).

Therefore, for **any** finite set of single-root fork tests, the response
matrix on the `k+1` classes `a = 0..k` has rank at most 2.  This is the exact
obstruction behind "distinct predictive classes do not imply linear
independence" in this family: the class count `k+1` is carried by a coordinate
`a` that the single-root language only ever sees linearly.

## 3. What raises the rank: independent roots (language depth), and nothing else

**Abstract model (exact).**  A bank of `k` gadgets, class `a` = number of
A-type gadgets, per-gadget clone-success fraction `h_A != h_B`.  A *d-root
test* picks `d` distinct gadgets (ordered, without replacement) and requires
each of the `d` roots/clones to succeed on its gadget.  The success
probability is

```
P_d(a) = sum_{t}  C(d,t) (a)_t (k-a)_{d-t} / (k)_d  * h_A^t h_B^{d-t},
```

a polynomial of exact degree `d` in `a` (verified by finite differences for
`k = 5`, `d = 1..5`).  The response matrix whose columns are the tests
`d = 0..k` (with `d = 0` the trivial test) has **exact rank `k+1`**
(Vandermonde-type), and the degeneracy control `h_A = h_B` collapses the rank
to 1 (the type coordinate is then unobservable at every depth).

**Answers to #599 Program B, at this level of resolution.**

* Q1 ("can `k+1` tests give rank `k+1` on the #549 classes?"): **within the
  single-root protocol, no** — rank ≤ 2 is exact (Proposition).  The abstract
  model shows the required language feature: `d >= k` **independent roots**
  (i.e. `d` sequential shared updates) give tests of degree up to `k` and a
  response matrix of rank `k+1`.
* Q3 ("smallest depth at which rank must grow?"): in this family **depth and
  rank growth are the same coordinate** — a language of depth `d` (d
  independent shared-update roots) has response rank at most `d+1` and can
  attain it; separating `k+1` classes needs depth `k`.  Clone *width* at fixed
  depth does not help (Proposition).
* Q2 ("can complete-survival language have small rank while branching grows?"):
  yes, and the mechanism is now explicit: the unbranched language is exactly
  the `d = 0` (depth-0) sector, rank ≤ 1 in the class coordinate, while
  depth-`d` branching languages reach rank `d+1`.

**Relation to #401's step structure.**  The same algebra appears there as the
`Z0/Z1/Z2/...` hierarchy: `q_m` needs `m`-point cooperative data and the exact
step-3 formula is cubic in `x` because it aggregates three-point events.  Root
count in the fork language plays the role of the interaction order in the
percolation language: response rank is bounded by the order of the language's
interaction algebra, and class count is not.

## 4. Boundary and the lift to the real N16 network

Everything above is exact but at the level of the *published protocol
constants*: per-gadget `EXIT` tables and the per-gadget success fractions of
the abstract model.  The abstract model assumes multi-gadget events factorize
per gadget.  The real #435 N16 network could satisfy or violate that
assumption at the level of *which* multi-root configurations are jointly safe.

**Cheapest falsification/verification of the lift** (no new sampling, one
exact enumeration): implement a two-root protocol on the actual #435 N16
pair — two sequential shared updates in *different* gadgets — and test whether
`P(a)` has a nonzero quadratic coefficient; the abstract model predicts it
does iff the joint-safe sets couple the two roots non-trivially.  This needs
the N16 configuration data (scripts of #435/#491), not new Monte Carlo.

**Claim boundary.**  Exact statements about the #549 protocol algebra and an
abstract exact gadget model; not a statement about the physical cut-network
threshold, not a claim that the abstract model equals the N16 network.

Files: `scripts/probe/c2_root_rank.py` (re-verifies #549 closed form and
produces the exact degree/rank table).
