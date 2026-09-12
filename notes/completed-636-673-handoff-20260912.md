# Updates for existing #636 / #673; no duplicate task

## #636 — the all-p width-four task is completed

Starting from #708's fixed input certificate, a 94-state common weighted
quotient and three explicit nonnegative local blocks of sizes 5,15,16 give
all-length polynomial identities:

```
P0 numerator = tr(B5^m)-tr(B16^m)+2t^(2m)
P2 numerator = tr(B15^m)-tr(B16^m)+2t^(2m)
M numerator  = tr(B15^m)-tr(B5^m)
```

Normalize by (1+t)^(4m), t=p/(1-p). The result is proved for all p, all m>=2
using 126 exact polynomial zero moments, a coefficient-bounded Kronecker
certificate and Cayley-Hamilton, not a parameter scan. Generic minimum scalar
orders are P0=17, P2=23, M=16; special p=1/2 gives 15, p=2/3 gives 14.

The leading Perron weights near the unique cylinder crossing are one on both
sides. The crossing has a degree-17 defining polynomial and agrees with the
already-published n=4 value q4=.59141717085313848.... A slower shared 0.773935
relative mode cancels identically from M; its visible relative decay is
0.251749755, with coefficient two. Rational certificates prove every finite
root lies below q4 and the fixed-width displacement formula in the note.

Do not re-dispatch this same visible-spectrum task or buy another width engine.
The genuinely larger question is an all-width closure/observable-weight theorem
and width-uniform bounds, not a list of more roots. No new hardware or scan is
commissioned by this comment. The thermal-derivative m*lambda^m control can be
cross-linked to #275 without changing that issue's candidate contract.

## #673 / #692 — rank-two onset and equality classification completed

The old shared-vertex union lower-bound argument is invalid, but its axis
conclusion is correct. Replace it with the two-cycle-core lemma. A minimal
ambient-rank-two core is a wedge or theta; a dumbbell is excluded by the
intersection form. Systolic lengths yield the axis onset 2L-1 and L^2 crosses.
For the diamond, all three theta cycle lengths are at least2L, yielding
3L-1 vertices; equality and period arithmetic force a full2L line with one
straight transverse L-1-site plug. There are exactly4L^2 minimizers for every
L>=2. Diamond L5 k14,count100 is a theorem, not a task needing2^50 enumeration.

Read `notes/two-cycle-core-onsets-20260912.md`; preserve the original note and
result history but remove the invalid proof dependency. No additional census
or new duplicate issue is needed. This does not classify every higher-mass
wrapping cell.

## Validation and integration boundary

Ten new local mathematical tests passed. The all-p proof, resultant, Hankel
minors and rational root inequalities were actually executed. A clean patch
application check is included in VALIDATION.json. Full repository CI has not
been run for this addition. No remote branch or issue was modified by the
analysis session. The parametric scripts use #708's existing certificate;
the onset proof/script is independent of #708 and may be integrated separately.
