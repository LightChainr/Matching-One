# One homogeneous N50 calculation of the original observer

This P1 calculation addresses the existing missing endpoint: epsilon=1,
t=0, on exactly the original parent pair (5,5)/(1,7), with periods
(a,b),(-b,a). All 50 sites are iid Bernoulli(p). Keep q=r-1, E=q^2 and
S*=2 beta-r-3K+2N+1=CB+CW+F4+Bvac in their existing units.
The source perturbation is exp(t S*), normalized separately in each
geometry. No new observable, descriptor, source fit, hole-order grid,
Monte Carlo block, m point or continuum-field attribution is introduced.

## Fixed question and interpretation

Compute the original common-root U and V=d_t U, including root and slope
motion, from the full finite population. The saturated parent endpoint
has V approximately +0.3891471785; its sign does not follow to epsilon=1
from the existing 0/1-hole tables. Record the fixed endpoint-sign
continuation prediction V>0 and the finite zero-transmission null V=0.
If the homogeneous rational interval is strictly negative, retire that
positive-sign continuation. If positive, reject the finite zero null
but do not call this a mechanism confirmation. An interval containing
zero is unresolved at the fixed arithmetic budget. Report U regardless
of its sign. There is no post-result extension to N100, different t,
epsilon or another source to repair a failed prediction.

This is deterministic finite-lattice calibration, not a new independent
statistical evidence block, an explanation of the asymptotic anomaly,
or reinstatement of the retired P154/P334/F4 production lines.

## Exact population and scoring

Use the full integer table (K,q,count,sum_S) for each geometry. At each K,
sum_q count must equal binomial(50,K), and the total must equal 2^50.
The sums retain qS and ES exactly because q is fixed within each row.
No conditional population is substituted for the complete law.

At t=0 use h=p/(1-p) and weights h^K. The common normalizer is (1+h)^50;
source derivatives still use each geometry's own S mean. Fix
DeltaCos4=-1152/625 and A50=50^(13/8)/2. Set M=(q1+q2)/2,
Y=(E1-E2)/DeltaCos4 and z=log(h). Evaluate at M=0:

    U/A50 = Y_z/M_z
    z_t = -M_t/M_z
    V/A50 = (Y_zt+z_t Y_zz)/M_z
            -Y_z*(M_zt+z_t M_zz)/M_z^2.

Use exact Fraction arithmetic with the already vendored interval class.
Require one Descartes sign variation of the pooled numerator, an opposite
sign bracket h in [0,4], exactly 160 rational bisections unless an exact
root is reached, and a strictly positive slope enclosure. All decisions
use rational intervals; A50-scaled decimals are display only. Failure of
a gate stops this score and is reported, with no extra precision search.

## Feasibility and stop

The lifted black-connectivity frontier has already reproduced the full
N9/N13 direct homology oracle and both existing N25 histograms. This does
not establish its full N50 resource budget. Initial N50 CPU/RSS probes
may stop before completion; they do not constitute a scientific score.
Compile a compact exact implementation and check its complete N25 tables
against the independently enumerated a70eeff0 tables before target use.
Freeze the implementation before any complete target run. Each N50 run
has a hard declared CPU/state/RSS gate, recorded with the output; a
partial run is not a histogram. If the measured peak remains uncontrolled,
publish the resource limit and do not start an open-ended cluster job.
No cloud startup is implied by this contract.

Commit this contract and score.py before completing a target table or
evaluating its root/U/V. Later implementation and resource-only changes
must preserve this contract and be identified in the run receipt.
