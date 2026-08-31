# Baseline reweighting is a measured part of the one-hole U response

Decision: **weighted_rank_jump_only_rejected** on the fixed N50 parent pair.

| Mixed original-U contribution | Exact coefficients evaluated numerically |
|---|---:|
| Baseline reweighting | +4.55032712323679 |
| Weighted observable jump | -15.3060455308009 |
| Total, imported from the completed one-hole result | -10.7557184075641 |

The two prescribed terms oppose each other: baseline reweighting partially
offsets the negative weighted jump. The one-term model misses a real positive
contribution; the full normalized defect operator already contains it, without
an added source or an adjustable mixing coefficient.

The reweighting contribution has a rational enclosure excluding zero:
True. This is the prescribed contribution
of Cov(w,O_intact), not a newly fitted residual or a rank-preserving population
share. Both rank-changing and rank-preserving configurations can contribute.
The jump contribution is total minus reweighting, so the two add exactly by
definition and do not count as independent evidence.

## How the missing information was obtained

The prior full defect packet lacks Sminus*qplus and Sminus*Eplus. The exact
single-defect topology restricts q/E changes to alternating four-neighbor
patterns. For every other configuration Oplus=Ominus; hence

`full sum(Sminus*Oplus) = old full sum(Sminus*Ominus)`
`+ alternating sum[Sminus*(Oplus-Ominus)]`.

Only 2 alternating patterns times 2^21 remaining bits were enumerated per
geometry: one eighth of the full population. Intact and defective observers
use identical B configurations. The Bernoulli degree stays25, including the
two forced occupied and two forced vacant neighbors. No baseline enumeration,
random sampling, root search or test suite was performed. Raw subset sums,
input commit/hash pointers and the exact arithmetic result are included.

## Fixed observer, source and chart

Source Sstar=C+F4+Bvac; pA=s+(1-s)p, pB=p, epsilon=1-s. At zero source the
baseline-reweighting insertion vanishes identically for every p. Its mixed
jet is hO=25(1-p)Cov(Sminus-Splus,Oplus). Each geometry is normalized before
pooling; the saved complementary root is reused. All four terms in
Xi_reweight/A = hY_p/D - Y_pp*hQ/D^2 - Y_p*hQ_p/D^2 + Y_p*Q_pp*hQ/D^3
are included, with D=Q_p and A=50^(13/8)/2. The derivative of25(1-p) is retained.

This decides the fixed weighted-rank-jump-only response model; it does not
identify a continuum field, finite-interior law or an independent production
effect. The earlier source-independent gain rejection, larger-N F4 unresolved
result and P154/P334 fixed decisions remain unchanged. The two contributions
are not two new adjustable sources. Do not fit another relative coefficient
to restore a failed one-term model.

Definition and mechanism algebra: notes/decimation-closed-source-and-global-u.md
at7132f0c2. Operator proof: bc17b81d:notes/checkerboard-single-defect-source.md.
Prior total Xi: f5c4a74a:results/p337-endpoint-defect/score/score.json.
Run: python scripts/analyze_defect_reweight.py --output-dir NEW_DIRECTORY.
