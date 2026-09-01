# Fixed N25, m64 two-law result

Freeze: `375a6f0ce67d46871ec97ea338fdf1342ed33e30`.

Primary decision: **certified_Ustar_negative_Udrop_positive_at_fixed_m64**.

The displayed bounds are outward-rounded rational enclosures; numerical midpoints do not decide signs. A25=25^(13/8)/2 is positive, so U and U/A25 have the same sign. Each law uses its own pooled root.

| Law | U/A25 lower | U/A25 upper | Sign |
| --- | --- | --- | --- |
| star | -62326039019195048123/10000000000000000000000000000000000000000 | -31163019509597524061/5000000000000000000000000000000000000000 | negative |
| drop | 11460271838603870300751659/10000000000000000000000000000000000000000 | 573013591930193515037583/500000000000000000000000000000000000000 | positive |

| Law / geometry | P(rank1) lower | P(rank1) upper | Necessary draws lower bound |
| --- | --- | --- | --- |
| star / axis | 5579946465849791227970987/5000000000000000000000000000000000000000 | 446395717267983298237679/400000000000000000000000000000000000000 | 851262647244163 |
| star / tilted | 106205815601482749261/5000000000000000000000000000000000000000 | 212411631202965498523/10000000000000000000000000000000000000000 | 44724481169877526072 |
| drop / axis | 20792294022790496528842449/5000000000000000000000000000000000000000 | 41584588045580993057684899/10000000000000000000000000000000000000000 | 228450020704475 |
| drop / tilted | 25123377123372850451/1000000000000000000000000000000000000000 | 251233771233728504511/10000000000000000000000000000000000000000 | 37813387719925333756 |

Resource gate: **do_not_promote_ordinary_unconditional_sampling_to_P0_new_estimator_required**.

The union bound uses ceil((19/20)/P1_upper). It is necessary for a 95% chance of seeing even one rank-one draw, regardless of dependence. It does not give a sufficient budget to estimate U or a lower bound for importance, conditional or other variance-reduced estimators. No wall-clock estimate follows.

The JSON includes both primitive integer root numerators, Descartes certificates, exact root brackets, positive slope enclosures, full-normalization observer derivatives, and all four population cells. No enumeration, simulation, cloud job, source fit or other coupling point was evaluated.

This calculation is a deterministic consequence of the existing exact histogram. It is not independent evidence, a continuum or large-N result, a uniform remainder bound, or a homogeneous continuation theorem. Failure of finite-point separation would not refute the proved eventual asymptotic signs.

Reproduce from the frozen checkout (standard-library Python >=3.10):

```sh
python3 experiments/p337-finite-law-window-20260831/score.py \
  --freeze-commit 375a6f0ce67d46871ec97ea338fdf1342ed33e30 \
  --output-dir /tmp/p337-m64-fresh-output
```
