# Fixed prospective safe-insertion estimands

One new square-site L512 insertion-clock block only: N=262144,
a=154646, b=155385, m=N-b=106759, q=4 probes per safe prefix,
B=14 independent equal-size batches of 10000 naturally sampled prefixes.
The scorer reads a completed `data/run.json` and its declared compressed
batches. This document and executable alone are not a production result.

## Prefixes, cohorts, and probes

Risk is `rank_b=1`. Early means `J1<=a`; late means `a<J1<=b`.
Let `c=nu_b`, `s=m-c`. For each safe prefix (`s>0`), each probe chooses
a uniformly random safe site in a fresh copy of that same prefix and records
`Y=nu(A+v)-c`. Average the four probes within prefix:

    Ybar_i = (Y_i0 + Y_i1 + Y_i2 + Y_i3)/4.

The probes are continuations of one prefix, not four independent prefixes.
Sampling must preserve the uniform-safe-site marginal even if probes repeat
a site. The input does not reveal probe locations or prove uniformity.
For `c=m`, no safe insertion exists; all four Y values must be -1 and Ybar
is undefined. These dead prefixes are excluded from primary and included
in secondary. Outside rank one, direction is (0,0), nu and all Y are -1;
rank zero has J1=0, rank two has 1<=J1<=b.

Rank-one D is the exact primitive unoriented integer direction: sign alone
is canonicalized so its first nonzero coordinate is positive. Nonprimitive
directions fail validation. Safe probes satisfy `0<=Y<=m-1-c`.

## Primary and its exact consumers

For each exact `(D,c)` cell h with s>0 and both cohorts present, use
prefix counts n_Eh,n_Lh and cohort means of prefix Ybar:

    w_h = n_Eh*n_Lh/(n_Eh+n_Lh),  W = sum_common w_h,
    dY_h = mean_Eh(Ybar) - mean_Lh(Ybar),
    deltaY = sum_common w_h*dY_h / W,
    weighted_earlyY = sum_common w_h*mean_Eh(Ybar) / W,
    weighted_lateY  = sum_common w_h*mean_Lh(Ybar) / W.

For the completion graph on the s safe vacancies, with e unordered edges,
the exact conditional mean is `E[Y|A]=2e/s`. Therefore within a matched cell:

    dE_h  = s_h*dY_h/2,
    dZ2_h = -s_h*dY_h/[m*(m-1)],
    deltaE  = sum_common w_h*dE_h/W,
    deltaZ2 = sum_common w_h*dZ2_h/W.

The last quantity estimates the early-minus-late **unconditional two-next-
insertion rank-one survival** contrast within matched cells, using
`P(J2>b+2|A)=[s(s-1)-2e]/[m(m-1)]`. It is not a directly observed binary
endpoint or survival conditional on the first step. Apply the s-dependent
factors inside cells before weighting. Both consumers use precisely the
primary support and weights. DeltaE and deltaZ2 are derived predictions
from the same probes, not additional independent measurements.

## Prespecified secondary replication

Use every rank-one prefix, including c=m. Within exact D only, form
`v_d=n_Ed*n_Ld/(n_Ed+n_Ld)` on shared cells and compute

    secondary_deltaNu = sum_common v_d*(mean_Ed(c)-mean_Ld(c))/sum_common v_d.

This repeats the old fixed direction-only completion-count contrast on a
new independent block. Its support and weights are separately recomputed;
it shares this new block with the primary. No old/new pooling occurs and
differences between these differently weighted targets are not mediated
fractions or explained percentages.

## Joint uncertainty, support, and outputs

The ordered vector is `(deltaY,deltaE,deltaZ2,secondary_deltaNu,
weighted_earlyY,weighted_lateY)`. Delete an entire aligned batch j from all
six coordinates, then recompute cell counts, support, means and weights:

    theta_bar = (1/B) sum_j theta_(-j),
    Cov_JK = (B-1)/B sum_j (theta_(-j)-theta_bar)(theta_(-j)-theta_bar)^T.

Report pooled plug-in point estimates, all 14 deletion vectors, full 6x6
covariance and its diagonal square-root SEs. The exact relations
`deltaY=weighted_earlyY-weighted_lateY` and
`deltaZ2=-2*deltaE/[m(m-1)]` make this covariance singular. No independence,
matrix inverse, or separate-probe standard error is assumed.

If either comparison has no common support in the full block or any
deletion, write `not_scoreable` diagnostics with all deletion vectors,
null the entire joint covariance/SEs, and exit 2. Undefined quantities
remain null; no deletion or cohort is silently dropped. Shared support may
otherwise change across deletions. Coverage uses both eligible prefixes
and all rank-one prefixes as denominators, separately by cohort and combined.
The overlap-weight sum is not asserted to be an effective sample size.

JSON retains cell counts, exact integer probe sums and squared sums,
prefix sums of squares, means, one-sided cells, support coverage, nonzero
probe/prefix counts, per-batch totals, input hashes/manifests, and scoring
sources. Production status, contract values, unique batch identities/files,
hashes, row counts, directions, and all sentinel/range rules are checked.
Structural errors exit 1 without scientific results. Input checks do not
validate the geometry engine or establish independence of generation.

Run after the parent has reviewed the scorer and acquired real data:

    python analyze.py
    python analyze.py --data-dir data --output-dir results-new
    python analyze.py --output results-new/result.json --md results-new/RESULT.md

Default inputs and outputs are relative to the script directory; an explicit
data directory defaults to `contract.json` in its parent. Explicit relative
CLI paths are relative to the invocation directory. `--contract` can name
the input contract. Default outputs are `results/result.json` and
`results/RESULT.md`; `--output-json`/`--output-md` are aliases. Existing
outputs are refused; select fresh paths for another scoring pass.

This primary examines one moment of one safe-step successor law. A
population discrepancy challenges recursive `(rank,D,nu)` sufficiency at
this fixed count; a weighted null may hide cell cancellation and does not
establish even full successor-law equality, much less recursive closure.
There is no subgroup search, cutoff search, fitted exponent, enlargement of
the acquired block, label-time conclusion, or continuum claim.
