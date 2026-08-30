# P250 radius-six parameter-free flat-extension freeze

## Question

The radius-five result fixed one cross-hand convention:
`Alexander R2 + coefficient conjugation`.  This follow-up does not reopen the
identity/R0/R1/R3 vote.  It asks the sharper question that radius five could not
answer: does the degree-two five-state quotient close as a genuinely flat
degree-three moment system, and does the selected map carry the whole
five-dimensional relation space rather than only its first null line?

For each hand form the two-charge Hankel block

```text
H3[(charge,u),v] = G_charge(u+v),   |u|,|v| <= 3,
```

with shape `20 x 10`.  The primary null is `rank(H3)<=5`.  If both hands pass,
compare the two five-dimensional right-kernel projectors after the already
selected R2-Alexander coordinate map and coefficient conjugation.

## Minimal sufficient acquisition

All entries with total degree at most five already exist in the independent
80k diamond and 1.2M radius-five archives.  The only missing block is
`degree 3 x degree 3`.  Its distinct source endpoints are

```text
(6,0), (5,1), (4,2), (3,3), (2,4), (1,5), (0,6).
```

For the minus hand, the fixed map is `phi(a,b)=(-a,b)`, so the target endpoints
are their quadrant-II images.  The spatial union has 13 points, but only 14
hand-point pairs: seven plus-source and seven minus-target.  With charges 1/2
this is 28 complex or 56 real coordinates per batch.

The full radius-six shell would acquire 96 complex coordinates when both hands
are retained.  The fixed-gauge subset is only `29.2%` of that payload.  Every
retained coordinate fills an entry of one of the two missing Hankel blocks;
the other 11 shell points fill neither.  This is the minimum sufficient subset,
not a convenience truncation.

The exact gate also confirms that all 13 displacements remain distinct on the
norm-101 parent quotient and avoid the origin.

## Frozen score

Maximum-volume pivot charts for ranks 5 through 9 are selected entirely from
the existing `12 x 10` top block, before any degree-six observation.  Their
condition numbers range from 88 to 325.  Each rank null uses its fixed Schur
complement `S-R P^-1 Q`; separate delete-one covariance contributions from the
80k, radius-five, and fresh streams are added.

The result logic is:

1. score `rank(H3)<=5` separately for plus and R2-gauged minus;
2. if either fails, use the frozen rank 6--9 ladder only to localize the
   truncated state-dimension lower bound;
3. only if both rank-five nulls survive, score the basis-independent kernel
   projector identity `conj(P_plus)=P_minus`;
4. never rescore the five old cross-hand candidates.

This is a path-independent commuting-moment test.  Endpoint moments still do
not distinguish ordered `TxTy` from `TyTx`, so a pass is not a magnetic or
noncommutative translation theorem.

## Budget

Freeze a fresh `1,200,000` replicas, 400 batches, 16 workers, new seed
`25060610120261250`, counters `[0,1200000)`.  Do not run from this commit.

The radius-five 1.2M job took 869 seconds on the same 16-core Huawei class.
The new payload is smaller, although common union-find and context construction
still dominate, so this remains a short server job.

At 1.2M the empirical radius-five coordinate standard-error proxy has median
`6.54e-5`.  The old-stream flat-completion predictor uncertainty is already
larger: median `2.12e-4` for plus and `1.53e-4` for minus.  Consequently the
fresh stream supplies only about 9% and 15% of median total variance.  Doubling
fresh samples would reduce total standard error by only about 2.2% and 4.0%.
The 1.2M choice reaches the inherited-covariance floor without wasting another
machine-hour.

For 28--56 resolved modes and alpha `0.01`, 80% power requires a noncentrality
of about 32--42, equivalent to an RMS whitened departure of `1.07--0.86` per
mode.  Thus this budget is aimed at an order-one predictor-standard-error
sixth-state departure, not at arbitrarily tiny violations.

## Scientific card

- Mechanism space: tests whether the selected radius-five line morphism lifts
  to a five-state flat quotient and a full relation-space bridge.
- Does not prove: ordered-path commutation, a microscopic quotient isomorphism,
  or a continuum operator identity.
- Observer/sector/source/geometry: fixed-p Z5 projective-leg two-charge moments;
  plus and minus Hecke hands; minimal norm-505 degree-six endpoint blocks.
- Dependency group: 80k degree-four archive, independent 1.2M radius-five
  archive, and one future independent radius-six stream.
- Natural lift: a pass freezes a five-state commuting quotient in the R2 gauge;
  a failure promotes the first surviving frozen rank as the truncated state
  lower bound, without reopening harmonic or morphism selection.
