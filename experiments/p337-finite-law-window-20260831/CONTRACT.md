# One fixed finite-coupling comparison and ordinary-sampling feasibility

Compare the already named laws Sstar=C+F4+Bvac and Sdrop=Sstar+r, with coefficient of r fixed at exactly one. Their opposite strong-coupling tails are proved, but no finite point for this two-law comparison has been scored here. The physical source and original q/E observer do not change. This is a deterministic consequence of an existing full finite histogram, not independent evidence.

## Fixed point and primary question

Use only N25 axis (5,0) and tilted (4,3), DeltaCos4=1152/625. Choose one multiplier m=64: it is the smallest dyadic m with N/m^2<=1/100. This scale rule motivates a finite strong-coupling probe; it is NOT itself a remainder bound or proof that the leading tail is accurate. The answer must come from the complete counts, not a truncated asymptotic expression.

There is no search for a peak, crossover or favorable multiplier. Do not append 64 to the older frozen [2,4,8,16] experiment. This is a separate two-law question with exactly one point. If it does not separate, report that outcome without increasing m or changing the rank coefficient.

With h=p/((1-p)*m), the exact weights are:

- star: count*h^K*m^(-g)
- drop: count*h^K*m^(-g+r), r=q+1.

For each law, normalize each geometry separately, find its own pooled root Q=(q_axis+q_tilt)/2=0, and compute the same original observer:

U/A25 = ((E_h,axis-E_h,tilt)/DeltaCos4) / ((q_h,axis+q_h,tilt)/2), A25=25^(13/8)/2.

The common p-to-h Jacobian cancels exactly. Do not compare the two laws at one law's root. Require an exact unique positive root (one Descartes sign variation of the pooled-root numerator), endpoint sign bracket and a strictly positive rational slope enclosure. Use 160 rational bisections, no adaptive precision or extra coupling grid. Verify the count at every K is binomial(25,K). Preserve complete normalization.

The primary decision is whether rational bounds strictly give Ustar<0<Udrop at this one finite coupling. Otherwise report the fixed-point comparison unresolved/inconsistent with that finite prediction, as applicable. Do not call absence of separation a disproof of the already established eventual asymptotic signs.

## Resource consequence tied to the original observer

Also compute each geometry's exact probability P1 of rank one, directly from its rank-one polynomial. This is required because E=1-P1; a raw unconditional sample with no rank-one configurations cannot directly resolve its thermal derivative.

For M ordinary draws from a specified law at its fixed root, the union bound gives Pr(any rank-one draw)<=M*P1, without independence assumptions. A necessary sample-count lower bound for a 95% chance of even one such draw is therefore ceil((19/20)/P1_upper). It is NOT a sufficient budget for estimating U, and does not bound conditional, importance, twist-partition or other variance-reduced estimators. No wall-clock forecast is inferred from sample count alone.

If this necessary bound exceeds 10^9 for any of the four law/geometry populations, do not promote direct unconditional sampling of this fixed comparison to P0. Retain the exact sign result and explicitly identify a new estimator as the missing requirement. No simulation is launched by this contract; 10 servers do not remove the rare-sector issue.

## Provenance and stop

All input Git blobs and hashes are in inputs/SOURCES.json. Commit this contract, input copies, dependency and score.py before the first root/U/P1 evaluation. No new enumeration, MC, cloud job, source fit, new observable or extension to another N or m. Report every law/geometry cell. The N25 Smith classes differ; no continuum, large-N uniformity or N50 saturation-to-homogeneous continuation follows from this fixed finite calculation.
