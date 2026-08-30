# Single-root theorem and projective ULC conjecture

Let `A_k` count size-`k` subsets whose homology image is rank one in a fixed
primitive line or stabilizer orbit, and normalize by the Boolean layer:

```text
q_k = A_k / binom(N,k).
```

The degree-`N-1` Bernstein coefficient of `A'(p)` is exactly

```text
d_k = (k+1) A_(k+1) - (N-k) A_k
    = (N-k) binom(N,k) (q_(k+1)-q_k).
```

Hence derivative signs are the discrete slope signs of the conditional layer
probability `q_k`.

## Exact minimal theorem

If the nonzero values of `q_(k+1)-q_k` change sign exactly once, from positive
to negative, then `A'(p)` has exactly one simple root in `(0,1)`. Strict
unimodality is not required.

Proof: set `t=p/(1-p)`. After removing the positive factor
`(1-p)^(N-1)`, the derivative is a polynomial `P_+(t)-P_-(t)` whose positive
coefficient degrees are all below its negative coefficient degrees. The ratio
`P_-/P_+` is strictly increasing: its logarithmic derivative is the difference
between degree averages on two disjoint ordered supports. It goes from zero to
infinity, crosses one exactly once, and has nonzero derivative there.

Contiguous support plus strictly decreasing adjacent ratios
`q_(k+1)/q_k` is a convenient sufficient condition. This is precisely strict
ULC of the raw layer counts on their positive support.

## Exact evidence and first qualification

All 16 focused orbit rows from the six-HNF atlas plus N13/N17 are strict single
peaks, strict ULC and strict ratio-decreasing.

The wider honest connected HNF scan through N12 contains 240 fixed-line and
217 stabilizer-orbit sequences. Every one has contiguous support, strict ULC,
strictly decreasing ratios and a single nonzero slope-sign change. But strict
single-peakedness holds only for 228/240 lines and 205/217 orbits. The first
plateau is N6 `[[2,0],[0,3]]`:

```text
q_2,q_3,q_4 = 1/5, 3/5, 3/5.
```

The plateau contributes one zero derivative coefficient; the exact theorem
still gives a unique simple root. Thus the right invariant is ULC/single
crossing, not a unique discrete mode.

## Why the obvious proof routes do not close

For a fixed line `ell`, the rank-one indicator is exactly

```text
1{H1=ell} = 1{ell subset H1} - 1{rank(H1)=2}.
```

Both right-hand terms are increasing events, and the fixed-line family is
order-convex because homology images grow under inclusion. But neither fact
alone proves ULC. A difference of nested upsets can have a normalized zero
valley; the four-element example `U\\V={{1},{2,3,4}}` already fails
log-concavity.

Nor is the site-ground homology rank a matroid rank. On the N4 quotient
`[[2,0],[0,2]]`, adding site 0 to sites `{1,2}` jumps ambient rank directly
from zero to two, violating the matroid unit-increment axiom.

## Research hierarchy

- **Theorem:** one normalized slope-sign change implies one simple balance root.
- **Exact bounded evidence:** strict ULC and ratio decrease for all 240 lines / 217 orbits through HNF N12, plus the focused N13/N17 rows.
- **Conjecture:** every honest finite square-torus fixed-line rank-one family has a contiguous ULC normalized layer sequence.
- **Proof target:** show these order-convex families have a normalized-matching or Lorentzian rank-generating property; a literal matroid or generic nested-upset theorem is insufficient.
