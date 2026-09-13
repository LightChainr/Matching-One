# Probe #624 re-run (2026-09-08) — location without shape: verdicts re-adjudicated on main

2026-09-08. Analysis + numpy/mpmath. No percolation engine, no arXiv, no
N=725, no #612 rescoring, no fit to #582 (1.55 retired, irrelevant here).
Does **not** enter `docs/STATUS.md`. Closes nothing. Issue #624 stays open.

## Provenance and what this PR is

The #624 probe was first delivered in **PR #634** (merged into
`claude/matching-one-workspace-pwr5pv`; files already present on main):

- `notes/probe-location-without-shape-20260907.md` — the theory note
- `scripts/probe/toy_cdf_families.py` — the toy families
- `results/probe-location-without-shape/latest.json` — the committed grid

The 2026-09-08 assignment re-opened the delivery channel with an override:
**draft PR against `main`**, not the old workspace branch. This PR is that
delivery. The mathematics is unchanged (nothing was wrong); what is new:

1. **Bit-exact reproduction on main.** `toy_cdf_families.py` rerun on main at
   this PR's head regenerates `latest.json` byte-identically (`git diff`
   empty). The committed artifact is reproducible on the target branch.
2. **Independent re-adjudication** (`toy_cdf_families_rerun_20260908.py`,
   new in this PR): every verdict re-derived through different code paths —
   mpmath at 50 dps, off the committed `u`-grid, at `N = 2^30` (beyond the
   committed `k ≤ 17`), and from the CDF side by inverting `F_N`. Six
   independent checks, all passing
   (`results/probe-location-without-shape/rerun-20260908.json`).
3. The committed `latest.json` records frontier
   `claude/matching-one-workspace-pwr5pv @ 8b5f9d1a` (its origin); this note
   is the provenance bridge to main.

## Verdicts (unchanged, now re-confirmed independently)

```text
T2  ACHIEVED (the prize). Z_N -> prescribed zeta; in these toys Z_N = zeta
    EXACTLY at every N. Re-confirmed at 50 dps at off-grid u and N = 2^30
    (i1: max dev 1.4e-42; i2 CDF-side: 4.1e-13, interp noise).
T1  ACHIEVED. Parity switch (gap 0.412, closed form i3) and
    t_N = (1+sin log N)/2 blend (t-spread 0.985 over the committed k).
T3  NULL. The anchor affine-orbit identity is EXACTLY zero at 50 dps
    (i5: 5.3e-51, i.e. rounding only — the float 0.0 in the original was
    not luck). Anchors are gauge for Route-A families.
D1  non-example re-confirmed. Route-B fixed-profile drift is 0.0 in the
    window for seven further decades of w_N (i4), matching exp(-c/w).
W5  numerical companion re-confirmed. Q-action invariance is exact for a
    second (alpha,beta) = (2.3,-0.77) (dev 3.9e-15); a p-axis kink at a
    different location (u_c = 0.35) and slope (2) still moves the shape
    O(1) (0.138) at every N.
```

## The one-line theory (for readers landing here first)

Master construction (Route A): for any continuous strictly increasing
`zeta` with `zeta(a)=0`, `zeta(b)=1` (anchors `a=0.2`, `b=0.8`), set

```text
Q_N(u) = p_* + w_N zeta(u),   w_N -> 0,   F_N := Q_N^{-1}.
```

Then `F_N` is a genuine CDF family with `F_N(0)=0`, `F_N(1)=1`, strictly
increasing about `p_*`, all interior quantiles in a `w_N`-window about `p_*`
(so `Q_N -> p_*` uniformly on every compact of `(0,1)` — Theorem L's
conclusion, cited from `notes/p613-quantile-convergence-20260907.md`,
not reproved), while `Z_N = zeta` **exactly** for every `N`. Location is
`w_N`'s business; shape is `zeta`'s; the two do not interact. Theorem L
cannot pin the shape: any conclusion of its strength is compatible with any
prescribed `Z_infinity` (T2), with no `Z_infinity` at all (T1), and with
anchors as pure gauge (T3-null).

Full statements and proofs: `notes/probe-location-without-shape-20260907.md`
(§2 construction + location proof, §3 T2 proposition, §4 T1, §5 D1,
§6 T3/D4 anchor gauge, §7 the two Aff(1) quotients, §8 D5 table). The D5
table (RSW / four-arm / self-duality / scaling-limit / finite-L: what pins
profile vs width) is unchanged and lives there; nothing in this re-run
touches its Y/N/unknown entries.

## Interface

- **#622 W1**: the analytic counterweight stands, now with a main-branch
  provenance chain and an independent verification artifact. Definition of
  percolation `S` remains #622's own job.
- **#618**: may cite Prop §3 of the 2026-09-07 note ("location does not
  imply a common `omega` for the 9-vector") — T2 is proved and machine-
  checked at 50 dps.
- **#613** cited (Theorem L), not reproved. **#620** not involved. No
  STATUS entry; issue left open; nothing merged or closed by this PR.

## Files (this PR)

```text
notes/probe-location-without-shape-20260908.md                    (this note)
scripts/probe/toy_cdf_families_rerun_20260908.py                  (independent re-adjudication)
results/probe-location-without-shape/rerun-20260908.json          (6/6 checks pass)
```

Run:

```text
/Users/lc/.workbuddy/binaries/python/envs/default/bin/python \
    scripts/probe/toy_cdf_families_rerun_20260908.py
```

(numpy + mpmath only; regenerates the rerun JSON. The original
`toy_cdf_families.py` rerun regenerates `latest.json` byte-identically.)
