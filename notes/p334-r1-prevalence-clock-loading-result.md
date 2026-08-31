# The integrated R1 contrast is almost a risk-state step

Final source `0d1e586dafbade5e7d1f9bfc598170d0c881e337` contains all
20,000 paired counters per N and 20 original batches. The prepared C/L
decomposition was applied once. Four-state covariance uses those same hybrid
`Y` records, with all 47 N325 and 164 N425 whole-pair fallbacks retained.
There is no new simulation or reliability solve.

## The main mechanism: between-risk-state variance

For `G=(Y_first-Y_second)/delta_cos4`, evaluated separately at `p_ref` and
after integration, condition on the four states `Rpair=(Rf,Rs)`:

\[
\operatorname{Var}(G)=
\operatorname{Var}(E[G\mid Rpair])+E[\operatorname{Var}(G\mid Rpair)].
\]

| N / endpoint | Between-state variance | Within-state variance | Between fraction ± batch SE |
|---|---:|---:|---:|
| 325 / canonical | 0.02651655 | 0.00496121 | 84.2390% ± 0.2844 percentage points |
| 425 / canonical | 0.01959476 | 0.00399282 | 83.0724% ± 0.3056 percentage points |
| 325 / integral | 0.11807528 | 0.00017862 | **99.84895% ± 0.00375 percentage points** |
| 425 / integral | 0.08802166 | 0.00012005 | **99.86380% ± 0.00374 percentage points** |

These are exact ddof-zero variance decompositions of the empirical per-counter
hybrid contrast. Their uncertainty comes from the same 20 original batches;
they are not variances of the overall estimated mean. N325 integral conditional
means in states `(00,01,10,11)` are approximately
`(0,+0.487407,-0.487185,+0.000557)`; N425 gives
`(0,+0.420679,-0.420936,-0.000214)`. Thus the large fluctuations are nearly
the steps between “only first at risk”, “only second at risk”, and the two
near-zero cases. This is not a claim that the clock carries no mean signal.

The exact integrated-clock identity explains the scale separation. With
`d=N-k0` and the declared hybrid waiting-time readout,

\[
m_i^{\rm int}=\frac{d+1-E[T_i\mid R1]}{N+1}.
\]

The observed conditional integrals are about `0.372` (N325) and `0.376`
(N425). R1 membership therefore introduces an order-0.37 step, while
within-R1 waiting-time variation is divided by `N+1`. The source's removed
suffix-noise fraction can accordingly be about 49–50% for canonical H4 but
under 1% for its integral. This does not contradict the previous 147-prefix
conditional-stratum 83.95% result: that conditioning, empirical mixture, and
readout variance denominator are different. It was never a global gain.

## Mean prevalence versus conditional loading

The original symmetric identity `D=C+L` remains separate from the variance
decomposition. Both terms below are H4 normalized, and their covariance is
retained:

| N / endpoint | C: prevalence ± SE | L: conditional clock ± SE | D=C+L |
|---|---:|---:|---:|
| 325 / canonical | −0.00035908 ± 0.00150262 | +0.00096906 ± 0.00039612 | +0.00060998 |
| 425 / canonical | +0.00168943 ± 0.00102780 | −0.00055355 ± 0.00043189 | +0.00113588 |
| 325 / integral | −0.00075551 ± 0.00316134 | +0.00016733 ± 0.00008074 | −0.00058818 |
| 425 / integral | +0.00359666 ± 0.00218668 | −0.00010548 ± 0.00007512 | +0.00349118 |

The signed mean terms partly cancel; prevalence-dominated **variance** does
not mean that the prevalence term dominates every **mean**. Equivalently,
`C_int=(rf-rs)[d+1-(ETf+ETs)/2]/[(N+1)delta]` and
`L_int=-(rf+rs)(ETf-ETs)/[2(N+1)delta]`. No causal or percentage attribution
is inferred.

The result saves the full 27-coordinate joint C/L and four-state LOO
covariance, rank at most 19 per size. No high-dimensional inverse covariance,
new omnibus test, or combination of correlated C/L readouts is used. This is
the R1 second-birth contribution only, not full `F2` or `A_top`. All results
reuse the original source blocks and keep their paired hybrid policy.

Scientific card: the key explanation is an almost discrete R1-membership
step behind integrated-contrast noise. The next microscopic decomposition
can separate original direct triggers, collective completion, and unchanged
fallback; it must remain in these same batches and the same covariance block.
Outputs: `results/p334-r1-prevalence-clock-loading/{score.json,REPORT.md}`.
