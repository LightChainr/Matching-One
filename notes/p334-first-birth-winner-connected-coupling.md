# Marginal winner allocation already reverses the marked source direction

All four integrated source point-direction reversals found at `2dd865f0`
remain after removing the connected K1/winner covariance. The connected term
instead weakens their magnitudes. The earlier reversal is therefore not
evidence that connected past/future coupling is its principal cause: the
accepted prevalence and conditional marginal first-birth/winner allocation
already produce it.

This is a descriptive decomposition of the same original populations.
Joint statistical uncertainty is delegated to the existing covariance
coordinator; none of these source signs is newly declared established.

## Exact decomposition on one fixed source population

For orientation i let I indicate **globally accepted and R1 in that
orientation**, exactly as in `2dd865f0`. All unconditional means below retain
the full20000-counter denominator. Put

```
r=E[I], k=E[I*K1], p=E[I*piD], m=E[I*mu],
H=E[I*H2*mu], B=E[I*K1*piD].
Kbar=k/r, pibar=p/r, c=B/r-Kbar*pibar.
```

Here piD is the exact conditional probability that the eventual second
birth is at an original-checkpoint singleton trigger, and mu is the full
conditional remaining wait. K1 and piD vary across accepted prefixes.
The connected quantity is the conditional empirical-population covariance
`c=Cov(K1,piD|I=1)`, without a Bessel factor. It also equals
`Cov(K1,1{eventual source is D}|I=1)` in the conditional suffix ensemble,
by the tower identity. It is not a claim of non-Markovian path memory.

The marked direct integral is

\[
 A_D^{\rm int}=\frac{H-B}{N+1}
 =\underbrace{\frac{H}{N+1}}_{\text{completion}}
 -\underbrace{\frac{kp/r}{N+1}}_{M_D:\ \text{marginal debt}}
 -\underbrace{\frac{B-kp/r}{N+1}}_{C_D:\ \text{connected debt}}.
\]

Since piG=1-piD, the collective counterparts are

\[
 F_{2,G}^{\rm int}=\frac{(d+1)r-m-H}{N+1},\quad
 M_G=\frac{k-kp/r}{N+1},\quad C_G=-C_D.
\]

Thus the connected term **only reallocates** marked loading between direct
and collective sources. It vanishes from their sum identically in the full
population and in every delete-one population; the untouched signed
remainder then restores the same global hybrid A. This cancellation does
not mean that each marked covariance is zero or uninformative.

## Actual conditional coordinates

| Size/orientation | Accepted-R1 r | Mean K1 | Mean piD | Cov(K1,piD) |
|---|---:|---:|---:|---:|
| N325 first | .32560 | 181.733876 | .661608884 | .197415930 |
| N325 second | .32645 | 181.526574 | .663382292 | .169422168 |
| N425 first | .32060 | 238.533843 | .650892013 | .153382826 |
| N425 second | .32195 | 238.334369 | .645652922 | .172944998 |

The accepted counts remain6512/6529 and6412/6439. Both conditional
covariances are positive at both sizes, but their paired orientation
contrast has opposite signs across sizes. This describes the archived
point coordinates, not an independent significance assessment.

## Which part produces the source reversals?

All entries below are the H4-normalized first-minus-second integral
contrasts. Full A equals completion minus marginal debt minus connected debt.

| Size/source | Completion | Marginal debt | Connected debt | A without connected | Full marked A |
|---|---:|---:|---:|---:|---:|
| N325 direct | +.0005547302 | +.0006529828 | -.0000360436 | **-.0000982526** | -.0000622091 |
| N325 collective | -.0001240242 | -.0003042295 | +.0000360436 | **+.0001802053** | +.0001441617 |
| N425 direct | -.0003571092 | -.0006157200 | +.0000171018 | **+.0002586109** | +.0002415090 |
| N425 collective | +.0008358653 | +.0012934718 | -.0000171018 | **-.0004576065** | -.0004405047 |

Every bold marginal-only entry already has the opposite sign from its
completion term. Restoring the connected term changes none of those signs
and moves each full marked contrast toward zero. Marginal-only here is an
algebraic product substitution within the fixed accepted population, not a
simulated independence intervention or a new fitted physical model.

To resolve the marginal term further while keeping prevalence separate,
write the unnormalized orientation difference as

\[
 \Delta(r\bar K\bar\pi)
 =\Delta r\,\overline{K\pi}
 +\bar r\,\bar\pi\,\Delta K
 +\bar r\,\bar K\,\Delta\pi.
\]

Bars in this display mean equal averages of the two orientation-conditional
coordinates; `overline(K*pi)` averages their two products. This fixed nested
midpoint identity first separates prevalence, then the two conditional
marginals. Divide by `(N+1)*delta_cos4` for the following direct-source debt:

| Size | Prevalence | Conditional K1 | Winner share | Marginal debt sum |
|---|---:|---:|---:|---:|
| N325 | +.0004109494 | -.0001799021 | +.0004219356 | +.0006529828 |
| N425 | +.0005485904 | -.0001092217 | -.0010550887 | -.0006157200 |

At N325, prevalence and winner-share terms reinforce, partly offset by the
K1 term. At N425, the winner-share term is large enough to outweigh the
opposing prevalence contribution. The conditional K1 term has the same
negative sign in both rows. This locates the main point-estimate change in
the orientation allocation of winning probability, not an enlarged connected
past/future covariance. The collective decomposition is saved alongside it.

## Sufficient-statistic and covariance handoff

Only six means per orientation are needed:
`[r,E(I*K1),E(I*piD),E(I*mu),E(I*H2*mu),E(I*K1*piD)]`.
The last four are read from the **already stored** source-clock,
positive-F2 and first-birth-subtraction batch vectors of `2dd865f0`.
Original gate states plus full births `9c495ab1` supply only r and E(I*K1).
No safe polynomial is evaluated, no permutation replayed and no DP/MC run.

For the point estimate, pool all20 sufficient-statistic batches first, then
calculate conditional ratios, products and covariances. For each LOO, drop
one original1000-counter batch, pool the remaining19 and recompute those
nonlinear quantities. Averaging per-batch products would define a different
quantity and is not used.

`results/p334-first-birth-winner-connected-coupling/score.json` retains the
20-by-12 `joint_20_batch_sufficient_means`, complete labels, pooled means,
the derived point vector and all20 correctly re-pooled LOO vectors. The
original source counts and full20k denominator are explicit. The derived
50-coordinate block is supplied for alignment with the global/source
covariance; no inverse or separate omnibus score is constructed here.

Scientific card: this separates a genuine marked past/future covariance
from conditional marginal allocation in the full-topology source signal.
It changes the explanation of the earlier point reversals while leaving
the complete A total and the fixed accepted population unchanged. Everything
shares the original e81dd59f stream, exact clocks0d1e586d and births9c495ab1;
there is no new independent evidence block or claim of population source
signs. The next useful readout is the common-covariance uncertainty on these
already named terms, not another model family or simulation.

```sh
/Users/lc/python-envs/research-py311/bin/python scripts/p334_first_birth_winner_connected_coupling.py
```
