# C10: ranked bridge conjectures and three executed cheap tests

**Status:** the probe's bridge-conjecture ranking, plus the three cheapest
tests executed in this round (T1 dynamic signatures on P398; T2 marked-channel
visibility of the C7 hidden coupling; T3 = re-verified protocol algebra, see
C2 — no new run needed).

## 1. Executed cheap tests

### T1 — does time evolution sharpen an even readout language to the orbit quotient?
Prediction (from round-2 Theorem 2 / C4 F1): an *even* language can never
separate two states in one C2 orbit, so its dynamic resolution is bounded by
`n_orbit`; the open question was how close it gets.  Measured dynamic
signature class counts on P398 (D0/D2 readouts, fingerprints over time
grids):

| w | states | C2 orbits | D0 static | D0 short-T | D0 mid/long-T | D2 short-T |
|---|---|---|---|---|---|---|
| 4 | 14 | 10 | 8 | 10 | 10 | 10 |
| 5 | 42 | 26 | 12 | 26 | 26 | 26 |
| 6 | 132 | 76 | 18 | 76 | 76 | 76 |
| 7 | 429 | 232 | 24 | 232 | 232 | 232 |
| 8 | 1430 | 750 | 32 | 746 | **750** | 750 |

(data: `results/probe-c10-t1-dynamic-signatures/latest.json`).
**Result:** time evolution sharpens the even dictionary to (essentially) the
orbit quotient at every width; the C9 unknown #3 ("does the D2 static gap
close dynamically?") is answered: **yes, up to the orbit ceiling** — and it
cannot exceed it.  The gap `209 -> 750` at w8 is closed by dynamics within a
short horizon.

### T2 — is the C7 hidden coupling visible through a marked (odd) channel?
Prediction: no even task sees `c`; a task with an odd source *and* an odd
readout should.  Measured (m=2,3; response at fixed times):

| channel | behaviour under c |
|---|---|
| even source + even readout | response columns identical for c = 0,.25,.5,1 (to 1e-15) — **invisible** |
| odd source + odd readout | response curves change strongly with c (c=0 frozen at its initial value; larger c decays faster) — **visible** |

(data: `results/probe-c10-t2-marked-channel/latest.json`).
**Result:** the marked-channel criterion of round 2 / #598 is confirmed on the
no-go family: visibility requires a nontrivial-character channel on the task
side; the hidden sector is not "weakly measurable", it is *exactly
invisible* to unmarked tasks.

### T3 — #549 protocol algebra (re-verified in C2; no new run)
The single-root fork probability is affine in `a` (rank ≤ 2), depth-d roots
give degree d (rank d+1 attainable).  This is the language-side counterpart of
T1/T2: resolution is bounded by what the declared language can express.

## 2. Bridge-conjecture ranking (with cheapest falsification)

| rank | conjecture | state | cheapest falsification / next action |
|---|---|---|---|
| 1 | Exact r_lin(w=8) exists and satisfies r_pos(w8) > r_lin(w8) | unknown (definition-blocked) | owner publishes #593's defining matrix; one mod-p elimination settles it |
| 2 | r_mem(w9) numerical order 15..16, effective order stays 4 (sublinear growth, energy saturates) | unknown | port frozen projected-memory construction + continuity gate; numpy run is cheap (C3 cost probe) |
| 3 | Even-language dynamic resolution = orbit quotient for every dictionary (T1 generalised) | **tested (holds w4–8, D0/D2)** | try a still-richer even dictionary or larger width to find the first gap |
| 4 | "One marked readout suffices to expose an odd perturbation" (naive reading of #598) | **false** (round 1: odd readout + even sources = invisible; T2 confirms channel pairs) | — |
| 5 | Two-root fork probability on the real N16 network has a nonzero quadratic coefficient in a | unknown | exact two-root enumeration on the #435 pair (needs N16 site data, no sampling) |
| 6 | p_c is visible to some current finite declared task | unknown (after C7, burden on a sector-completeness bridge) | check whether declared readouts see the homology/topological sector that E_p[X]=0 targets (#581/#584 typed labels) |
| 7 | A smooth one-parameter family can mimic all #584 remainder labels | unknown | run the discriminator suite (cocycle/curvature/crossed labels) on archived Q_N data (needs data access) |
| 8 | Compositional depth ≤ d ⇔ algebraic degree ≤ d for general tree languages | unknown | formal definition of compositional depth; then test on #550's depth-2 witness |

## 3. What the ranking says

The top four items are now either answered (4 = false, 3 = true on tested
dictionaries) or blocked on definitions/data (1, 6, 7) or on a mechanical
port (2).  Item 5 is the best "open science" target that needs no new
sampling and no owner input beyond the N16 site data already in the repo.

Related: C9 ladder (`probe-complexity-ladder-v1-20260906.md`); C2
(`probe-fork-root-count-algebra`); C7 (`probe-threshold-no-go-theorem`).
