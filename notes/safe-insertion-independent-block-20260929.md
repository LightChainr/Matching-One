# Independent L512 safe-insertion block: unresolved drift, unreplicated direction signal

2026-09-29. A new, fixed 140,000-prefix square-site L512 block, **not** a
reanalysis of the preceding 140k block. Acquisition source and contract were
pushed at `55fa265629701e04b879b774ae7d3b97fefa845f` before sampling. The
scorer was pushed at `f340aa55134213be71d22941fae2ccb7fd068cb7` while the
job was running and before any production response was read.

**Outcome.** The prespecified completion-creation contrast is unresolved.
The prespecified independent replication of the earlier direction-only
completion-count contrast also does **not** reproduce that signal. This
lowers the empirical standing of the earlier retrospective direction
challenge; it does not establish a Markov model or invalidate the exact
finite L4 counterexamples.

## 1. Actual new acquisition and primary readout

Keep the old insertion counts `a=154646`, `b=155385`, N=262144. Among
rank-one prefixes at b, early means `J1<=a`, late means `a<J1<=b`.
For each prefix, four independent with-replacement choices from its safe
vacancies are inserted in **separate copies**. Their completion-count
increments Y are averaged within prefix. A rejected immediately completing
vacancy is resampled; this gives the uniform safe-site law. Probe random
numbers use a separate stream from the site permutation.

The primary contrast compares prefix-average Y at the same exact primitive
direction D and integer completion count nu. Cell overlap weights are
`n_E*n_L/(n_E+n_L)`. The secondary compares current nu at the same D, under
its own weights, repeating the earlier chosen statistic independently.

| Prespecified quantity, early minus late | Estimate | Aligned batch SE |
|---|---:|---:|
| Primary: safe-insertion completion increment, same (D,nu) | +0.00668183 | 0.01208706 |
| Implied mean synergy-edge contrast, same (D,nu) | +356.3978 | 645.0460 |
| Implied two-next-insertion survival contrast, same (D,nu) | −6.25403e−8 | 1.13192e−7 |
| Secondary: current completion count, same D | −0.0398651 | 0.5698779 |

Weighted early/late means of Y are 0.0956368 and 0.0889549. All errors are
14-whole-batch delete-one SEs, with support and weights recomputed after
each deletion. The full six-dimensional covariance and all deletion vectors
are saved. Four probes are **not** counted as four independent prefixes.

The middle two rows are exact-map estimators from the same probes, not
independent outcomes or observed binary exits. For m=N-b, s=m-nu in each
cell, the [pair-synergy identity](completion-pair-synergy-20260929.md) gives

    Delta(mean e) = s*Delta(mean Y)/2,
    Delta(two-step survival) = -s*Delta(mean Y)/[m(m-1)].

Apply s inside each cell before averaging. None of these two-step quantities
replaces the old 735-insertion survival endpoint. Different weighting targets
do not define an explained/mediated fraction.

## 2. The independent replication changes the scientific interpretation

The earlier [retrospective within-direction readout](directional-hazard-contrast-20260929.md)
was **+1.62917 +/- 0.39098** completion sites. The same defined comparison
on this new block is **−0.03987 +/- 0.56988**. Old and new blocks are not
pooled. The prior positive point estimate cannot continue to be presented
as a stable empirical exclusion of direction-only closure.

Both observations remain in the record. The new result neither proves zero
nor explains why the old selected comparison was larger. In particular,
there is no basis here to assert a discovered engine error, a physical
sign reversal, or a universal disappearance of memory. Its proper status
is **retrospective signal not independently reproduced**.

For the primary, the new block has not resolved a difference in the mean
creation of completion sites after matching D and nu. Even a precise zero
of that weighted mean would not establish the entire successor law, exclude
opposing cell effects, or prove recursive closure. These finite count-clock
observations do not settle the continuum law.

The exact L4 [weak-Markov counterexamples](completion-pair-synergy-L4-20260929.md)
remain true. Their existence disproves universal finite-size closure; it
does not guarantee a detectable history effect at L512 or under this one
weighted moment. Likewise the earlier L3 closure remains valid at L3.

## 3. What the acquisition actually resolves

Of 140,000 prefixes, 55,028 are rank one at b: 24,201 early and 30,827
late. No rank-one prefix has all remaining sites immediately completing.
The primary's 1,190 common cells retain 54,463 prefixes, **98.973%** of
rank-one risk. Missing overlap is reported, not filled by binning or a
different conditional definition.

There are 220,112 valid probes but only 1,233 nonzero increments (about
0.560%); 1,215 prefixes have at least one nonzero probe. Thus the simulator
has measured many independent prefixes but most continuation measurements
return zero. This sparsity is a directly observed limitation of this probe
design; it alone does not quantify a variance decomposition or prove the
entire SE is continuation noise.

The experiment therefore redirects attention:

- Treat L512 direction-only closure as unresolved, with an unreplicated
  retrospective challenge. Do not keep using the old point estimate as an
  accepted constraint on the next mechanism.
- The next numerical gain should come from a **variance-reduced measurement
  of the same completion-creation target**, for example analytic conditional
  averaging or an unbiased, explicitly weighted sampler. An unweighted
  convenience sample of high-contact sites changes the source and is not a
  replacement. No such sampler is claimed implemented here.
- Another full small-size census would not settle the large-size issue;
  merely extending these sparse uniform probes has less decision value than
  reducing their noise or deriving a continuum-kernel prediction. Theory on
  the joint near-critical law remains the parallel primary consumer.

This is attention allocation, not a prohibition on parallel ideas. No extra
samples, time cuts, directions, feature fits or response-selected subgroup
analyses were added after seeing the result.

## 4. Execution and reproduction

Only TV2N0X (`4a8d1d443419434889e49148ed0a7ba6`) was operated, in the new
directory `/workspace/Matching-One-TV2N0X/safe-insertion-20260929`.
Its 14.5-CPU / 25-GiB actual cgroup budget was measured, and 14 workers used.
The restarted ARM64 container required GCC installation; GCC 10.3.1 and
Python 3.9.9 ran the acquisition. No other machine, data directory or old
production block was altered.

The 32-prefix excluded benchmark forecast 353.99 s of ideal-parallel work;
actual complete acquisition was **519.843 s** (14 x 10,000). The benchmark
screen is an estimate, not a promised wall time. The declared sample count
was retained; there was no response-dependent extension. A 20-site probe
control on the two known physical L4 witness configurations ran once, not
a new census. The scorer had a small in-memory arithmetic check, not a
full test suite; the main agent read it before its one complete data pass.

All raw compressed rows, seeds, source/binary/batch hashes and timestamps
are in [data/run.json](../analysis/safe-insertion-production-20260929/data/run.json).
[Result JSON](../analysis/safe-insertion-production-20260929/results/result.json)
retains full covariance and per-cell sufficient statistics;
[RESULT.md](../analysis/safe-insertion-production-20260929/results/RESULT.md)
is the compact numerical report. Results were scored locally in the existing
research Python 3.11 environment, with no new local dependencies.

```bash
python3 analysis/safe-insertion-production-20260929/analyze.py \
  --output-dir /tmp/matching-safe-score-new
```

The scorer preserves existing outputs; choose an unused output directory.
Do not rerun the acquisition simply to reproduce the score. The remote
raw files are retained, while the completed job's machine is returned to
Ready after confirming no other research processes remain.
