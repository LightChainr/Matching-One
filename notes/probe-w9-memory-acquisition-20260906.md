# C3: w9 memory-order computation — cost probe and acquisition decision

**Status:** cost probe done; the full w9 projected-memory run is recommended
as the next bounded heavy step (not run here to keep this session within
probe discipline).  No science beyond the budget is claimed.

## Measured cost (numpy sparse uniformization; see `scripts/probe/c3_cost_probe.py`)

| width | states | nnz | rate | Poisson terms @ t=2 | one e^{tG}v @ t=2 |
|---|---|---|---|---|---|
| 8  | 1 430 | 14 014 | 56 | 203 | 0.014 s |
| 9  | 4 862 | 53 768 | 64 | 225 | 0.059 s |
| 10 | 16 796 | 206 856 | 72 | 246 | 0.244 s |

State enumeration and generator build are negligible (≤0.4 s at w10).

## Consequences

1. **The repository's "widths 9–10 are unreachable in pure Python" is a code
   limit, not a hard barrier.**  With a sparse numpy stack a single matrix
   exponential at w9 costs ~0.06 s; even a few dozen such operations (the
   projected-memory/MZ kernel object set of #588) fit comfortably in a normal
   session.  w10 is ~4x more expensive per vector but remains feasible for
   vector-count workloads.
2. **What actually gates the w9 run is not compute but convention porting:**
   the frozen #580/#588 projected-memory construction
   (`scripts/p398_projected_memory.py`) must be ported off the repository's
   width-≤8 `Generator` cap and re-gated against the committed w8 numbers
   (the same continuity-gate pattern #603 used).  That port is a mechanical
   but real task and is the recommended next step if the repository wants the
   w9 saturation number before #593's C++ path lands.
3. **Decision table** (for whoever runs it):
   * numerical order (tol 1e-6) at w9 = 15..16 and effective 99.9% order
     stays 4 ⇒ "numerical order grows sublinearly, energy order saturates"
     (current reading extrapolated);
   * numerical order at w9 stays 14 ⇒ saturation confirmed earlier than
     expected;
   * effective order rises ⇒ #588's "BOUNDED_MEMORY_ORDER_IS_ACCURACY_
     DEPENDENT" needs the accuracy ladder re-examined.
4. The same stack also makes an exact w9 `r_lin`-style Krylov computation
   feasible once the #593 definition of `r_linear` is published (see
   `probe-rlinear-definition-calibration-20260906.md`).

## Boundary

Cost numbers are machine-specific but the *scaling* (nnz ~ 11–12 per row,
Poisson terms ~ +10% per width step, wall time per expv ~ x4 per width step)
is robust.  No memory-order scientific value is claimed from this probe
alone.
