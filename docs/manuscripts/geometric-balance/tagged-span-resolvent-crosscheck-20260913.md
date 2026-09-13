# Cross-check of the tagged-span resolvent delivery

2026-09-13. Verification note for the additive package handed over as
`matching_one_tagged_span_resolvent_20260913`, which is
`docs/manuscripts/geometric-balance/tagged-span-resolvent.md` plus
`span-resolvent-frontier.md` and four data/code files, and which was delivered
in the same round as comment 5653178247 on #739. Nothing in the delivered files
is changed here. This note records what we were able to check independently, and
one reproduction caveat that the delivered `EXECUTION.json` states too strongly.

Reproduce with

    python scripts/tagged_span_crosscheck.py --write

which reads the delivered `tagged-span-resolvent.json`, our own
`span-spectrum-20260913.json`, and `winding-nu-certified-20260913.json` (the
output of `scripts/winding_nu_certified.py`, an exact-rational winding-cluster
engine that shares no code with the delivered construction), and writes
`tagged-span-crosscheck-20260913.json`.

## 1. Reception of the package

All six delivered files reproduce `EXECUTION.json`'s sha256 exactly after
`git apply`, and the two machine-independent test runs agree with the handover:

| check | result |
|---|---|
| sha256 of all six files vs `EXECUTION.json` | 6/6 match |
| `python -m unittest discover -s tests -p 'test_tagged_winding_span.py'` | 17 tests, OK, 0.46 s |
| `scripts/tagged_span_controls.py --output ...` | 132994 bytes, 29.9 s class, all 14 systems closed |

Closing widths are 2--8 for both adjacencies, with tagged state counts
5, 13, 43, 131, 411, 1275, 3963 and exit-law lumps NN 3, 5, 10, 17, 36, 71, 161 /
matching 2, 3, 7, 15, 33, 68, 152, matching the table in section 4 of the note.
`build()` raises rather than truncating above width 8, so nothing here rests on
a silently incomplete automaton.

**Caveat.** `EXECUTION.json` records
`"report_regenerated_identically_except_elapsed_seconds": true`. That is true on
the delivering machine but not across platforms. Re-running the controls report
on macOS/py3.13.12 differs from the delivered file in two ways:

1. `laplace_mean_scaled_float_diagnostic` differs in the last one or two bits
   (relative ~1e-16). This is explicitly labelled a floating diagnostic and is
   immaterial.
2. For the widths 5--8 systems, the `centres` and all `*_interval` bounds are
   *different exact rationals with different denominators*, because a centre is
   a float solve plus one correction, and the float solve depends on the BLAS.
   Relative difference between the two platforms' centres: at most 9.34e-32.

That second item is not a defect, and it is worth saying why: the two platforms'
certified intervals overlap pairwise in 16/16 w>=5 systems, and — the real test
— **both platforms' intervals contain our independently computed exact rational
`nu_w` in all 12 systems where an independent exact value exists** (widths
6, 7, 8). A certificate that survives a different BLAS and an independently
written exact engine is doing its job. But "regenerated identically" should not
be quoted as a portability claim; the correct claim is "regenerated within its
own certified intervals".

## 2. Density certificate vs an independent exact value

`scripts/cylinder_winding_intensity.py` computes the winding-cluster density
`nu_w` by exact rational Gauss-Jordan on a transfer built from scratch, then
re-certifies the forward error through the empty-row reset. It shares no code
with `tagged_winding_span.py`. Comparing its exact value against the delivered
section 4.1 intervals:

| graph | w | p | interval contains our exact nu | interval width | centre rel. error |
|---|---:|---:|---:|---:|---:|
| NN | 6 | 1/8 | yes | 1.16e-31 | 7.5e-33 |
| NN | 6 | 1/4 | yes | 3.67e-31 | 3.4e-32 |
| NN | 7 | 1/8 | yes | 3.70e-31 | 4.2e-33 |
| NN | 7 | 1/4 | yes | 1.02e-30 | 2.1e-32 |
| NN | 8 | 1/8 | yes | 4.70e-31 | 5.2e-32 |
| NN | 8 | 1/4 | yes | 4.07e-30 | 8.4e-33 |
| matching | 6 | 1/16 | yes | 4.14e-32 | 2.2e-32 |
| matching | 6 | 1/8 | yes | 1.68e-31 | 3.4e-32 |
| matching | 7 | 1/16 | yes | 1.04e-31 | 9.5e-34 |
| matching | 7 | 1/8 | yes | 1.41e-31 | 4.2e-32 |
| matching | 8 | 1/16 | yes | 6.46e-31 | 1.9e-32 |
| matching | 8 | 1/8 | yes | 6.59e-31 | 2.5e-32 |

12/12 contained, no violations. The delivered centres sit ~1e-32 relative from
the exact value, which is the expected accuracy of one Newton-like correction on
a float solve; the intervals are a few orders looser than that but still
astronomically tight. Two independent exact engines now agree on `nu_w` at
w=6,7,8 for both adjacencies and five p values. The parameter-6 and
parameter-8 controls of the delivered `winding-prefactor` line are therefore
consistent, not merely self-consistent.

## 3. Our own depth-clamped spectrum against the all-height law

Our `span-spectrum-20260913.json` measures a depth-clamped chain with a D_MAX
bin cutoff plus a tail bin, so its `E[L]` and `CV^2` are censored. Comparing all
30 cells against the delivered all-height moments:

| tail fraction of cell | cells | worst CV^2 relative bias |
|---|---:|---:|
| <= 1e-6 | 26 | 2.4e-05 |
| > 1e-6 | 4 | 1.8e-01 |

The 4 materially censored cells are all at p=1/2, where the winding cluster is
not small on a narrow cylinder:

| graph | w | p | D_MAX | tail fraction | CV^2 truncated | CV^2 all-height | rel. bias |
|---|---:|---:|---:|---:|---:|---:|---:|
| NN | 4 | 1/2 | 48 | 1.57e-06 | 0.304422 | 0.304476 | 1.8e-04 |
| matching | 2 | 1/2 | 48 | 1.29e-06 | 0.559489 | 0.559600 | 2.0e-04 |
| matching | 3 | 1/2 | 48 | 2.11e-03 | 0.591669 | 0.625015 | 5.3e-02 |
| matching | 4 | 1/2 | 48 | 2.41e-02 | 0.523973 | 0.637515 | 1.8e-01 |

Consequences for what we already published. The erratum's section 2 already
flagged censoring in this file and added a per-run `censored` flag; this
quantifies it. Every cell used in the section 5.1 moment-limit revision is
either uncensored to <= 2.4e-05 relative, or — for the w=7 rows — carries a tail
fraction of order 1e-7 to 1e-12. The two p=1/2 rows that are off by 5% and 18%
in CV^2 were not used in the slope test, and must not be quoted as measurements
of the p=1/2 span law.

## 4. The two constructions agree height by height, not just in their moments

The strongest check available, and the one that actually used the delivered
construction rather than treating it as one more set of numbers. The two engines
compute the same object through structurally unrelated state spaces:

| | state space at w=7 | stored quantity |
|---|---:|---|
| `span_spectrum_build.cpp` | 389391 states (668439 at w=6) | all component ages, then projected onto a span histogram |
| `tagged_winding_span.py` | 71 lumped states (36 at w=6) | one tagged lineage, no age, no depth cutoff |

Agreement at every height h <= D_MAX and on the cutoff tail, across all 30 cells
of the committed spectrum:

| quantity | worst relative difference, over scored heights |
|---|---:|
| `d_h` at w=2,3,4 | 1.0e-14 |
| `d_h` at w=5 | 1.4e-12 |
| `d_h` at w=6 | 2.4e-12 |
| `d_h` at w=7 | 3.0e-12 |
| cutoff tail bin, w=5..7 | 3.3e-12 |

Two further details make this a real cross-validation rather than a coincidence.

First, the error grows smoothly with height and with width (4e-16 at h=11 to
3e-12 at h=48) — exactly the signature of the committed file's own float64
stationary solve, whose accuracy degrades as `d_h` becomes a small remainder.
The tagged side is exact rational, so the residual is the older file's error,
not a constructional disagreement. At w=2, both are exact and the difference is
identically zero.

Second, 231 of the 1408 compared heights fall below `1e-20 * nu`, where the
committed file stores denormal-scale floats (`d_h` ~ 1e-40, tail bins ~ 1e-45);
those are excluded from the score rather than used to characterise agreement, and
are reported separately as `heights_below_floor`. Without that floor the naive
worst case reads 0.73 purely from noise in the 1e-40 range. The check is
scale-aware for that reason.

A caveat on what this does and does not establish. Both engines could share a
misreading of what "span" means — the check cannot detect that, and the
arithmetic agreement would survive it. What it does establish is that the
depth-clamped construction and the tagged construction are computing the same
quantities to the accuracy of the less precise of the two, so the tag/forbidden
lumping is not silently changing the observable. A 1e-12-level agreement between
a 389391-state and a 71-state description of the same law is also the clearest
statement available of how much the age bookkeeping was buying.

## 5. The section 6.3 moment limit, now tested on uncensored moments

The point of the delivered construction for us is section 2's
`nu = alpha Z b`, `m1 = alpha Z^2 b`, `m2 = alpha (2Z^3 - Z^2) b`, which give
`E[L]` and `CV^2` at every closed width with no height cutoff. That extends the
moment-limit test from four families on truncated w=4..7 data to six families on
all-height w=2..8 data.

Write `T = pi/3 - 1 = 0.04719755119659775` and
`R_w = (CV^2 - T) * w`. Model A (`CV^2 = T + c/w`) makes `R_w` a constant.
Model B (`CV^2 = c/w`, i.e. the ratio tends to zero) forces `R_w` to fall at
slope exactly `-T` per unit width.

| graph | p | widths | R_w at w=4,5,6,7,8 | slope of R_w | rms(A) | rms(B) | verdict |
|---|---|---:|---|---:|---:|---:|---|
| NN | 1/8 | 4-8 | 0.66857, 0.72404, 0.74127, 0.73550, 0.71725 | +0.0109 | 0.0258 | 0.0847 | model A |
| NN | 1/4 | 4-8 | 0.65687, 0.66484, 0.65916, 0.64985, 0.63973 | -0.0049 | 0.0086 | 0.0600 | model A |
| matching | 1/16 | 4-8 | 0.27734, 0.28788, 0.29296, 0.28981, 0.28399 | +0.0015 | 0.0054 | 0.0691 | model A |
| matching | 1/8 | 4-8 | 0.52247, 0.54774, 0.55598, 0.55555, 0.55046 | +0.0064 | 0.0124 | 0.0762 | model A |

Model B requires `R_8 = R_4 - 4T = R_4 - 0.1888`, i.e. about 0.4798 / 0.4681 /
0.0886 / 0.3337. Measured: 0.7173 / 0.6397 / 0.2840 / 0.5505. The gap is
0.17--0.24 in `R_w` units, roughly an order of magnitude larger than the whole
observed drift of `R_w` across four widths, and model A's rms is 3--13 times
smaller than model B's in every family. Model B is excluded; model A is
consistent.

The two p=1/2 families have only widths 2--4 in the delivery, so they cannot
enter a slope test; they are excluded rather than reported as support.

Two things the table does *not* say. First, `R_w` is not flat: NN p=1/8 rises
from 0.6686 to 0.7413 by w=6 and then falls, and NN p=1/4 falls monotonically
after w=5. A one-parameter `c/w` correction cannot produce a turnover, so the
finite-width drift is not a clean single power and no correction exponent should
be quoted from it. Second, the gap at w=8 is between model A and model B, not
evidence for `T` itself: a family can fit model A well while `c(p)` is still far
from its limit. Section 6.3's conjecture remains neither proved nor refuted —
it is refuted only in the strong "ratio tends to zero" reading, which was
already the reading withdrawn in the erratum.

## 6. What changes, and what is now open

Changed: the w=8 point of the moment-limit test exists, exactly and without
censoring, so the w=8 span-table build is no longer on the critical path; it is
retained only as a third exact route to the same two numbers. The p=1/2 rows of
our committed spectrum are now known to be censored by 5% and 18% in CV^2 and
must be labelled as such wherever they are quoted.

Unchanged: the delivered note's section 5 criticisms of the read commit are
accurate. Items 3 and 4 (the `light=True` branch not dividing by `p_den**w`, and
the certificate being computed for the rounded candidate `a/2^k` rather than the
printed `pi`) are the same two defects our erratum found and fixed independently;
the delivered file `tagged-span-resolvent.json` was indeed absent from the read
commit's tree, which is the gap the erratum restored. Item 1 — that a complete
component's span is a maximum-minus-minimum, not a sum of irreducible piece
heights — targets a mechanism reading, not the computation; the erratum had
already withdrawn that mechanism, and the height-by-height agreement in section 4
above — where a 389391-state age-tracking chain and a 71-state tagged chain land on
the same `d_h` to 1e-12 — confirms the truncated chain was measuring the true span
all along.

Open, in rough order of value:

1. Width 9, past the delivered `build()` cap of 8, to see whether `R_w` keeps
   drifting or plateaus — the cleanest remaining discriminator among
   correction exponents.
2. A method independent of linear algebra. Every check above, including the
   height-by-height one, is exact rational or float64 arithmetic on one or
   another transfer operator; a direct simulation of the cylinder is the only
   check that does not share that architecture.
3. The width-8 height-by-height row of section 4, which is currently missing
   because no committed truncated w=8 spectrum exists. The delivered tagged
   construction can supply the reference side immediately — at NN p=1/4,
   `nu = 4.432636548604856e-05`, `E[L] = 4.56334396378038`,
   `CV^2 = 0.12716350006804397`, `sum_{h>20} d_h = 1.2246438200402218e-11`,
   and at matching p=1/8, `nu = 3.791735919415724e-05`,
   `E[L] = 4.711805413158092`, `CV^2 = 0.11600493483091046`,
   `sum_{h>20} d_h = 1.3400726036376886e-11`; both close exactly,
   `sum(d_h) + tail = nu`.
4. The `psi_w(z,u,v)` joint activity of the delivered section 3, whose
   `Corr(L,K)^2 ~ 0.867` at NN w=4 is the frontier note's conjecture J input.
   The delivered direct-activity transfer closes only to width 5 or 6, so this
   is where a width extension would buy something the resolvent does not
   already give.

Boundary. This note certifies arithmetic, mutual consistency, and platform
robustness of the delivered finite construction. It does not verify the
unique-anchor pathwise theorem, the completeness of the tagged automaton beyond
the recorded widths, or any asymptotic statement. "Model A consistent" means
consistent with the finite data, not established.
