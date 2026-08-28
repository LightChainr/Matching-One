# Issue #8 — L^{-4} / L^{-6} cancellation in matched estimators

Post-processing of Issue #7 exact torus microcanonical totals and Issue #9
Newman–Ziff matched roots. **No Monte Carlo was run.** Weights for
Richardson cancellation depend only on lattice sizes. Training residuals
are fit diagnostics and are **not** uncertainties. Intercept min–max across
models is **model spread**, not a statistical confidence interval.
Estimator ranking uses withheld-size prediction error only, never closeness
to a published threshold.

## Inputs

Issue #7 exact microcanonical files:

- `L02_microcanonical.csv`: L=2, occupancy k=0..4
- `L03_microcanonical.csv`: L=3, occupancy k=0..9
- `L04_microcanonical.csv`: L=4, occupancy k=0..16
- `L05_microcanonical.csv`: L=5, occupancy k=0..25

Issue #9 shared roots on L = 16, 24, 32, 48, 64, 96, 128, 192, 256
(independent baseline on L = 32, 64, 128, 256). Observables: wrapping
H, V, either, both, and cluster matching M. Shared-mode wrapping roots
coincide to printing precision because of the Issue #7 wrapping identity.
Cluster matching M is noisier; some large-L M roots left `[0.590, 0.595]`
and are used as numeric roots, not discarded.

## Protocol

Rolling one-step folds (shared mode):

- train through 64 → predict 96
- train through 96 → predict 128
- train through 128 → predict 192
- train through 192 → predict 256

On each fold the ordinary correction **model and L_min are chosen only on
training sizes** (inner one-step holdout of the last training size; ties
break toward the simpler model, then larger L_min). The chosen pair is
frozen, refit on all training sizes, then scored on the withheld size.

L^{-4} cancellation uses adjacent sizes with `w1+w2=1` and
`w1 L1^{-4}+w2 L2^{-4}=0`. L^{-4}/L^{-6} cancellation uses three consecutive
sizes and also annihilates L^{-6}. Combined uncertainty is the SD/SE of
the 20 batch combinations `p_super[b]=sum_i w_i p_{L_i}[b]`.

Matching annihilator (when microcanonical statistics exist):

```
A(p) = L1^{13/4} F_{L1}(p) - L2^{13/4} F_{L2}(p) = 0
```

with `F = M_L` for cluster matching and `F = D_L^x` for wrapping class x.
The annihilator root is scored against the ordinary withheld-size root.
Per-batch annihilator SD is not available: Issue #9 stored pooled
microcanonical sufficient statistics, not per-batch histograms.

Marks BETTER / SAME / WORSE compare **median absolute withheld error**
to the ordinary rolling estimator on the same observable. They describe
out-of-sample numbers only.

## Per-observable withheld-size comparison (shared, rolling-aligned)

Ordinary = nested model/L_min selection predicting `p_{L_next}`.
L^{-4}-cancel = last adjacent training pair, scored vs `p_{L_next}`.
L^{-4}/L^{-6}-cancel = last consecutive training triple, scored vs `p_{L_next}`.
Matching annihilator = `A(p)=0` on last adjacent training pair, scored vs ordinary `p_{L_next}`.

### Observable `H`

| estimator | n_success | median_abs_err | worst_abs_err | ratio vs ordinary | mark | noise amp (emp median / L1 median) |
|---|---:|---:|---:|---:|---|---|
| ordinary rolling | 4 | 5.7448217e-5 | 0.00014101803 | 1 | — | n/a |
| L^-4 cancel | 4 | 0.00015004649 | 0.00021662551 | 2.61186 | WORSE | emp 1.34418 / L1 1.70901 |
| L^-4/L^-6 cancel | 4 | 0.00017348703 | 0.00024182899 | 3.01989 | WORSE | emp 1.60616 / L1 2.19785 |
| matching annihilator | 4 | 0.00015054648 | 0.00021639175 | 2.62056 | WORSE | L1 2.01198 (no batch SD) |

Fold-level ordinary rolling:

| train_max | withheld_L | model | L_min | predicted | true | signed err | abs err | std err (err/SE) |
|---:|---:|---|---:|---|---|---|---|---|
| 64 | 96 | A | 32 | 0.592731906286681628 | 0.592872924315839933 | -0.000141018029158304640 | 0.000141018029158304640 | -0.894289937575844716 |
| 96 | 128 | A | 16 | 0.592782369166817075 | 0.592692551895070086 | 0.0000898172717469889845 | 0.0000898172717469889845 | 0.570388610815850895 |
| 128 | 192 | A | 24 | 0.592735842664208249 | 0.592760921825630005 | -0.0000250791614217562457 | 0.0000250791614217562457 | -0.123218147403258394 |
| 192 | 256 | A | 32 | 0.592761268006348679 | 0.592777028564663633 | -0.0000157605583149540902 | 0.0000157605583149540902 | -0.0785305414247944439 |

Fold-level L^{-4} cancel (last training pair):

| L1 | L2 | withheld_L | combined | true | abs err | emp noise amp | L1 amp |
|---:|---:|---:|---|---|---|---|---|
| 48 | 64 | 96 | 0.592724687926152793 | 0.592872924315839933 | 0.000148236389687139966 | 1.42507230173380146 | 1.92571428571428571 |
| 64 | 96 | 128 | 0.592909177408585604 | 0.592692551895070086 | 0.000216625513515517646 | 1.15725948355528304 | 1.49230769230769231 |
| 96 | 128 | 192 | 0.592609065231742328 | 0.592760921825630005 | 0.000151856593887676754 | 1.66673596165450003 | 1.92571428571428571 |
| 128 | 192 | 256 | 0.592777751346998600 | 0.592777028564663633 | 7.22782334967446154e-7 | 1.26329080249444318 | 1.49230769230769231 |

Fold-level L^{-4}/L^{-6} cancel (last training triple):

| L1 | L2 | L3 | withheld_L | combined | true | abs err | emp noise amp | sum\|w\| | cond |
|---:|---:|---:|---:|---|---|---|---|---|---|
| 32 | 48 | 64 | 96 | 0.592722651189544575 | 0.592872924315839933 | 0.000150273126295357951 | 1.74492768971578493 | 2.43645320197044335 | 35374949850.2 |
| 48 | 64 | 96 | 128 | 0.592934380889792272 | 0.592692551895070086 | 0.000241828994722185763 | 1.25600249853887779 | 1.95925058548009368 | 257443771103. |
| 64 | 96 | 128 | 192 | 0.592564220883478391 | 0.592760921825630005 | 0.000196700942151614458 | 2.01973061050967958 | 2.43645320197044335 | 2.26399560959e+12 |
| 96 | 128 | 192 | 256 | 0.592800795898263119 | 0.592777028564663633 | 0.0000237673335994855063 | 1.46739453819305806 | 1.95925058548009368 | 1.64763996481e+13 |

Fold-level matching annihilator:

| L1 | L2 | withheld_L | annihilator root | ordinary withheld | abs err | status | sum\|w\| |
|---:|---:|---:|---|---|---|---|---|
| 48 | 64 | 96 | 0.592724683895926141 | 0.592872924315839933 | 0.000148240419913791980 | ok | 2.29271725405110148 |
| 64 | 96 | 128 | 0.592908943644342878 | 0.592692551895070086 | 0.000216391749272792443 | ok | 1.73124766219742035 |
| 96 | 128 | 192 | 0.592608069291012729 | 0.592760921825630005 | 0.000152852534617276004 | ok | 2.29271725405110148 |
| 128 | 192 | 256 | 0.592777873444189085 | 0.592777028564663633 | 8.44879525452105904e-7 | ok | 1.73124766219742035 |

### Observable `V`

| estimator | n_success | median_abs_err | worst_abs_err | ratio vs ordinary | mark | noise amp (emp median / L1 median) |
|---|---:|---:|---:|---:|---|---|
| ordinary rolling | 4 | 5.7448217e-5 | 0.00014101803 | 1 | — | n/a |
| L^-4 cancel | 4 | 0.00015004649 | 0.00021662551 | 2.61186 | WORSE | emp 1.34418 / L1 1.70901 |
| L^-4/L^-6 cancel | 4 | 0.00017348703 | 0.00024182899 | 3.01989 | WORSE | emp 1.60616 / L1 2.19785 |
| matching annihilator | 4 | 0.00015054648 | 0.00021639175 | 2.62056 | WORSE | L1 2.01198 (no batch SD) |

Fold-level ordinary rolling:

| train_max | withheld_L | model | L_min | predicted | true | signed err | abs err | std err (err/SE) |
|---:|---:|---|---:|---|---|---|---|---|
| 64 | 96 | A | 32 | 0.592731906286681628 | 0.592872924315839933 | -0.000141018029158304640 | 0.000141018029158304640 | -0.894289937575843020 |
| 96 | 128 | A | 16 | 0.592782369166817075 | 0.592692551895070086 | 0.0000898172717469889845 | 0.0000898172717469889845 | 0.570388610815834501 |
| 128 | 192 | A | 24 | 0.592735842664208249 | 0.592760921825630005 | -0.0000250791614217562457 | 0.0000250791614217562457 | -0.123218147403256081 |
| 192 | 256 | A | 32 | 0.592761268006348679 | 0.592777028564663633 | -0.0000157605583149540902 | 0.0000157605583149540902 | -0.0785305414247958013 |

Fold-level L^{-4} cancel (last training pair):

| L1 | L2 | withheld_L | combined | true | abs err | emp noise amp | L1 amp |
|---:|---:|---:|---|---|---|---|---|
| 48 | 64 | 96 | 0.592724687926152793 | 0.592872924315839933 | 0.000148236389687139966 | 1.42507230173383956 | 1.92571428571428571 |
| 64 | 96 | 128 | 0.592909177408585604 | 0.592692551895070086 | 0.000216625513515517646 | 1.15725948355527182 | 1.49230769230769231 |
| 96 | 128 | 192 | 0.592609065231742328 | 0.592760921825630005 | 0.000151856593887676754 | 1.66673596165450370 | 1.92571428571428571 |
| 128 | 192 | 256 | 0.592777751346998600 | 0.592777028564663633 | 7.22782334967446154e-7 | 1.26329080249443832 | 1.49230769230769231 |

Fold-level L^{-4}/L^{-6} cancel (last training triple):

| L1 | L2 | L3 | withheld_L | combined | true | abs err | emp noise amp | sum\|w\| | cond |
|---:|---:|---:|---:|---|---|---|---|---|---|
| 32 | 48 | 64 | 96 | 0.592722651189544575 | 0.592872924315839933 | 0.000150273126295357951 | 1.74492768971584916 | 2.43645320197044335 | 35374949850.2 |
| 48 | 64 | 96 | 128 | 0.592934380889792272 | 0.592692551895070086 | 0.000241828994722185763 | 1.25600249853885436 | 1.95925058548009368 | 257443771103. |
| 64 | 96 | 128 | 192 | 0.592564220883478391 | 0.592760921825630005 | 0.000196700942151614458 | 2.01973061050968166 | 2.43645320197044335 | 2.26399560959e+12 |
| 96 | 128 | 192 | 256 | 0.592800795898263119 | 0.592777028564663633 | 0.0000237673335994855063 | 1.46739453819304809 | 1.95925058548009368 | 1.64763996481e+13 |

Fold-level matching annihilator:

| L1 | L2 | withheld_L | annihilator root | ordinary withheld | abs err | status | sum\|w\| |
|---:|---:|---:|---|---|---|---|---|
| 48 | 64 | 96 | 0.592724683895926141 | 0.592872924315839933 | 0.000148240419913791980 | ok | 2.29271725405110148 |
| 64 | 96 | 128 | 0.592908943644342878 | 0.592692551895070086 | 0.000216391749272792443 | ok | 1.73124766219742035 |
| 96 | 128 | 192 | 0.592608069291012729 | 0.592760921825630005 | 0.000152852534617276004 | ok | 2.29271725405110148 |
| 128 | 192 | 256 | 0.592777873444189085 | 0.592777028564663633 | 8.44879525452105904e-7 | ok | 1.73124766219742035 |

### Observable `either`

| estimator | n_success | median_abs_err | worst_abs_err | ratio vs ordinary | mark | noise amp (emp median / L1 median) |
|---|---:|---:|---:|---:|---|---|
| ordinary rolling | 4 | 5.7448217e-5 | 0.00014101803 | 1 | — | n/a |
| L^-4 cancel | 4 | 0.00015004649 | 0.00021662551 | 2.61186 | WORSE | emp 1.34418 / L1 1.70901 |
| L^-4/L^-6 cancel | 4 | 0.00017348703 | 0.00024182899 | 3.01989 | WORSE | emp 1.60616 / L1 2.19785 |
| matching annihilator | 4 | 0.00015054648 | 0.00021639175 | 2.62056 | WORSE | L1 2.01198 (no batch SD) |

Fold-level ordinary rolling:

| train_max | withheld_L | model | L_min | predicted | true | signed err | abs err | std err (err/SE) |
|---:|---:|---|---:|---|---|---|---|---|
| 64 | 96 | A | 32 | 0.592731906286681628 | 0.592872924315839933 | -0.000141018029158304640 | 0.000141018029158304640 | -0.894289937575833340 |
| 96 | 128 | A | 16 | 0.592782369166817075 | 0.592692551895070086 | 0.0000898172717469889845 | 0.0000898172717469889845 | 0.570388610815873085 |
| 128 | 192 | A | 24 | 0.592735842664208249 | 0.592760921825630005 | -0.0000250791614217562457 | 0.0000250791614217562457 | -0.123218147403255719 |
| 192 | 256 | A | 32 | 0.592761268006348679 | 0.592777028564663633 | -0.0000157605583149540902 | 0.0000157605583149540902 | -0.0785305414247951334 |

Fold-level L^{-4} cancel (last training pair):

| L1 | L2 | withheld_L | combined | true | abs err | emp noise amp | L1 amp |
|---:|---:|---:|---|---|---|---|---|
| 48 | 64 | 96 | 0.592724687926152793 | 0.592872924315839933 | 0.000148236389687139966 | 1.42507230173384175 | 1.92571428571428571 |
| 64 | 96 | 128 | 0.592909177408585604 | 0.592692551895070086 | 0.000216625513515517646 | 1.15725948355528755 | 1.49230769230769231 |
| 96 | 128 | 192 | 0.592609065231742328 | 0.592760921825630005 | 0.000151856593887676754 | 1.66673596165452875 | 1.92571428571428571 |
| 128 | 192 | 256 | 0.592777751346998600 | 0.592777028564663633 | 7.22782334967446154e-7 | 1.26329080249443264 | 1.49230769230769231 |

Fold-level L^{-4}/L^{-6} cancel (last training triple):

| L1 | L2 | L3 | withheld_L | combined | true | abs err | emp noise amp | sum\|w\| | cond |
|---:|---:|---:|---:|---|---|---|---|---|---|
| 32 | 48 | 64 | 96 | 0.592722651189544575 | 0.592872924315839933 | 0.000150273126295357951 | 1.74492768971584480 | 2.43645320197044335 | 35374949850.2 |
| 48 | 64 | 96 | 128 | 0.592934380889792272 | 0.592692551895070086 | 0.000241828994722185763 | 1.25600249853888507 | 1.95925058548009368 | 257443771103. |
| 64 | 96 | 128 | 192 | 0.592564220883478391 | 0.592760921825630005 | 0.000196700942151614458 | 2.01973061050972572 | 2.43645320197044335 | 2.26399560959e+12 |
| 96 | 128 | 192 | 256 | 0.592800795898263119 | 0.592777028564663633 | 0.0000237673335994855063 | 1.46739453819303465 | 1.95925058548009368 | 1.64763996481e+13 |

Fold-level matching annihilator:

| L1 | L2 | withheld_L | annihilator root | ordinary withheld | abs err | status | sum\|w\| |
|---:|---:|---:|---|---|---|---|---|
| 48 | 64 | 96 | 0.592724683895926141 | 0.592872924315839933 | 0.000148240419913791980 | ok | 2.29271725405110148 |
| 64 | 96 | 128 | 0.592908943644342878 | 0.592692551895070086 | 0.000216391749272792443 | ok | 1.73124766219742035 |
| 96 | 128 | 192 | 0.592608069291012951 | 0.592760921825630005 | 0.000152852534617053959 | ok | 2.29271725405110148 |
| 128 | 192 | 256 | 0.592777873444189085 | 0.592777028564663633 | 8.44879525452105904e-7 | ok | 1.73124766219742035 |

### Observable `both`

| estimator | n_success | median_abs_err | worst_abs_err | ratio vs ordinary | mark | noise amp (emp median / L1 median) |
|---|---:|---:|---:|---:|---|---|
| ordinary rolling | 4 | 5.7448217e-5 | 0.00014101803 | 1 | — | n/a |
| L^-4 cancel | 4 | 0.00015004649 | 0.00021662551 | 2.61186 | WORSE | emp 1.34418 / L1 1.70901 |
| L^-4/L^-6 cancel | 4 | 0.00017348703 | 0.00024182899 | 3.01989 | WORSE | emp 1.60616 / L1 2.19785 |
| matching annihilator | 4 | 0.00015054648 | 0.00021639175 | 2.62056 | WORSE | L1 2.01198 (no batch SD) |

Fold-level ordinary rolling:

| train_max | withheld_L | model | L_min | predicted | true | signed err | abs err | std err (err/SE) |
|---:|---:|---|---:|---|---|---|---|---|
| 64 | 96 | A | 32 | 0.592731906286681628 | 0.592872924315839933 | -0.000141018029158304640 | 0.000141018029158304640 | -0.894289937575843485 |
| 96 | 128 | A | 16 | 0.592782369166817075 | 0.592692551895070086 | 0.0000898172717469889845 | 0.0000898172717469889845 | 0.570388610815834501 |
| 128 | 192 | A | 24 | 0.592735842664208249 | 0.592760921825630005 | -0.0000250791614217562457 | 0.0000250791614217562457 | -0.123218147403256081 |
| 192 | 256 | A | 32 | 0.592761268006348679 | 0.592777028564663633 | -0.0000157605583149540902 | 0.0000157605583149540902 | -0.0785305414247947198 |

Fold-level L^{-4} cancel (last training pair):

| L1 | L2 | withheld_L | combined | true | abs err | emp noise amp | L1 amp |
|---:|---:|---:|---|---|---|---|---|
| 48 | 64 | 96 | 0.592724687926152793 | 0.592872924315839933 | 0.000148236389687139966 | 1.42507230173387226 | 1.92571428571428571 |
| 64 | 96 | 128 | 0.592909177408585604 | 0.592692551895070086 | 0.000216625513515517646 | 1.15725948355527134 | 1.49230769230769231 |
| 96 | 128 | 192 | 0.592609065231742328 | 0.592760921825630005 | 0.000151856593887676754 | 1.66673596165450043 | 1.92571428571428571 |
| 128 | 192 | 256 | 0.592777751346998600 | 0.592777028564663633 | 7.22782334967446154e-7 | 1.26329080249443832 | 1.49230769230769231 |

Fold-level L^{-4}/L^{-6} cancel (last training triple):

| L1 | L2 | L3 | withheld_L | combined | true | abs err | emp noise amp | sum\|w\| | cond |
|---:|---:|---:|---:|---|---|---|---|---|---|
| 32 | 48 | 64 | 96 | 0.592722651189544575 | 0.592872924315839933 | 0.000150273126295357951 | 1.74492768971589102 | 2.43645320197044335 | 35374949850.2 |
| 48 | 64 | 96 | 128 | 0.592934380889792272 | 0.592692551895070086 | 0.000241828994722185763 | 1.25600249853885314 | 1.95925058548009368 | 257443771103. |
| 64 | 96 | 128 | 192 | 0.592564220883478391 | 0.592760921825630005 | 0.000196700942151614458 | 2.01973061050967657 | 2.43645320197044335 | 2.26399560959e+12 |
| 96 | 128 | 192 | 256 | 0.592800795898263119 | 0.592777028564663633 | 0.0000237673335994855063 | 1.46739453819304789 | 1.95925058548009368 | 1.64763996481e+13 |

Fold-level matching annihilator:

| L1 | L2 | withheld_L | annihilator root | ordinary withheld | abs err | status | sum\|w\| |
|---:|---:|---:|---|---|---|---|---|
| 48 | 64 | 96 | 0.592724683895926141 | 0.592872924315839933 | 0.000148240419913791980 | ok | 2.29271725405110148 |
| 64 | 96 | 128 | 0.592908943644342878 | 0.592692551895070086 | 0.000216391749272792443 | ok | 1.73124766219742035 |
| 96 | 128 | 192 | 0.592608069291012951 | 0.592760921825630005 | 0.000152852534617053959 | ok | 2.29271725405110148 |
| 128 | 192 | 256 | 0.592777873444189085 | 0.592777028564663633 | 8.44879525452105904e-7 | ok | 1.73124766219742035 |

### Observable `M`

| estimator | n_success | median_abs_err | worst_abs_err | ratio vs ordinary | mark | noise amp (emp median / L1 median) |
|---|---:|---:|---:|---:|---|---|
| ordinary rolling | 4 | 0.0054819384 | 0.010233088 | 1 | — | n/a |
| L^-4 cancel | 4 | 0.0080023081 | 0.013443175 | 1.45976 | WORSE | emp 1.31298 / L1 1.70901 |
| L^-4/L^-6 cancel | 4 | 0.0088756121 | 0.014496046 | 1.61906 | WORSE | emp 1.63144 / L1 2.19785 |
| matching annihilator | 4 | 0.0078664556 | 0.016854608 | 1.43498 | WORSE | L1 2.01198 (no batch SD) |

Fold-level ordinary rolling:

| train_max | withheld_L | model | L_min | predicted | true | signed err | abs err | std err (err/SE) |
|---:|---:|---|---:|---|---|---|---|---|
| 64 | 96 | B | 16 | 0.594657884763017812 | 0.591867838280574210 | 0.00279004648244360179 | 0.00279004648244360179 | 1.26682935771469563 |
| 96 | 128 | A | 16 | 0.593417394926817474 | 0.586579496171288572 | 0.00683789875552890244 | 0.00683789875552890244 | 2.50451484874863675 |
| 128 | 192 | C | 32 | 0.587341836003921940 | 0.597574924071262892 | -0.0102330880673409517 | 0.0102330880673409517 | -3.02985756339447109 |
| 192 | 256 | A | 16 | 0.593001244711296340 | 0.588875266755668703 | 0.00412597795562763657 | 0.00412597795562763657 | 1.09698300420676231 |

Fold-level L^{-4} cancel (last training pair):

| L1 | L2 | withheld_L | combined | true | abs err | emp noise amp | L1 amp |
|---:|---:|---:|---|---|---|---|---|
| 48 | 64 | 96 | 0.594542459174626568 | 0.591867838280574210 | 0.00267462089405235764 | 1.48830198683159457 | 1.92571428571428571 |
| 64 | 96 | 128 | 0.591177888112514216 | 0.586579496171288572 | 0.00459839194122564391 | 1.13766523383045141 | 1.49230769230769231 |
| 96 | 128 | 192 | 0.584131749252133505 | 0.597574924071262892 | 0.0134431748191293867 | nan | 1.92571428571428571 |
| 128 | 192 | 256 | 0.600281490938948878 | 0.588875266755668703 | 0.0114062241832801755 | nan | 1.49230769230769231 |

Fold-level L^{-4}/L^{-6} cancel (last training triple):

| L1 | L2 | L3 | withheld_L | combined | true | abs err | emp noise amp | sum\|w\| | cond |
|---:|---:|---:|---:|---|---|---|---|---|---|
| 32 | 48 | 64 | 96 | 0.594402114759321722 | 0.591867838280574210 | 0.00253427647874751180 | 1.78653763013796314 | 2.43645320197044335 | 35374949850.2 |
| 48 | 64 | 96 | 128 | 0.590718247257034386 | 0.586579496171288572 | 0.00413875108574581444 | 1.47635164567129969 | 1.95925058548009368 | 257443771103. |
| 64 | 96 | 128 | 192 | 0.583078877928168571 | 0.597574924071262892 | 0.0144960461430943205 | nan | 2.43645320197044335 | 2.26399560959e+12 |
| 96 | 128 | 192 | 256 | 0.602487739803267919 | 0.588875266755668703 | 0.0136124730475992155 | nan | 1.95925058548009368 | 1.64763996481e+13 |

Fold-level matching annihilator:

| L1 | L2 | withheld_L | annihilator root | ordinary withheld | abs err | status | sum\|w\| |
|---:|---:|---:|---|---|---|---|---|
| 48 | 64 | 96 | 0.594545426340734817 | 0.591867838280574210 | 0.00267758806016060722 | ok | 2.29271725405110148 |
| 64 | 96 | 128 | 0.591182288035070114 | 0.586579496171288572 | 0.00460279186378154171 | ok | 1.73124766219742035 |
| 96 | 128 | 192 | 0.580720315746335203 | 0.597574924071262892 | 0.0168546083249276892 | outside_requested_interval | 2.29271725405110148 |
| 128 | 192 | 256 | 0.600005386127223872 | 0.588875266755668703 | 0.0111301193715551690 | outside_requested_interval | 1.73124766219742035 |

## All adjacent-pair / triple constructions (shared), not just rolling-aligned

These rows use every adjacent pair or consecutive triple and score against
the next larger size. They are extra diagnostics; the marks above use only
the rolling-aligned subset.

| observable | method | n_success | median_abs_err | worst_abs_err | median emp noise amp | median L1 amp |
|---|---|---:|---:|---:|---:|---:|
| H | L4 all pairs | 7 | 0.00014823639 | 0.00026849805 | 1.27178 | 1.70901 |
| H | L4L6 all triples | 6 | 0.00017348703 | 0.00024182899 | 1.46739 | 2.43645 |
| H | annihilator all pairs | 7 | 0.00014824042 | 0.00026834345 | n/a | 2.01198 |
| V | L4 all pairs | 7 | 0.00014823639 | 0.00026849805 | 1.27178 | 1.70901 |
| V | L4L6 all triples | 6 | 0.00017348703 | 0.00024182899 | 1.46739 | 2.43645 |
| V | annihilator all pairs | 7 | 0.00014824042 | 0.00026834345 | n/a | 2.01198 |
| either | L4 all pairs | 7 | 0.00014823639 | 0.00026849805 | 1.27178 | 1.70901 |
| either | L4L6 all triples | 6 | 0.00017348703 | 0.00024182899 | 1.46739 | 2.43645 |
| either | annihilator all pairs | 7 | 0.00014824042 | 0.00026834345 | n/a | 2.01198 |
| both | L4 all pairs | 7 | 0.00014823639 | 0.00026849805 | 1.27178 | 1.70901 |
| both | L4L6 all triples | 6 | 0.00017348703 | 0.00024182899 | 1.46739 | 2.43645 |
| both | annihilator all pairs | 7 | 0.00014824042 | 0.00026834345 | n/a | 2.01198 |
| M | L4 all pairs | 7 | 0.0026746209 | 0.013443175 | 1.38522 | 1.70901 |
| M | L4L6 all triples | 6 | 0.0033365138 | 0.014496046 | 1.63144 | 2.43645 |
| M | annihilator all pairs | 7 | 0.0026775881 | 0.016854608 | n/a | 2.01198 |

## Ordinary full-window intercepts (model spread, not a CI)

Fits on all shared sizes with L_min = 16. Training RMSE is a residual,
not an uncertainty.

| observable | model | pc | a4 | a6 | a8 | cond | train RMSE | withheld L=none |
|---|---|---|---|---|---|---|---|---|
| H | A | 0.592766703219904303 | -2.93089043260586597 | nan | nan | 211786.563238 | 0.0000729771273296119009 | full data |
| H | B | 0.592746545538598116 | 66.2275800717690729 | -17710.7763060976454 | nan | 478256245.242 | 0.0000628830069835224090 | full data |
| H | C | 0.592770923234743795 | -381.483870930549311 | 363823.590954631032 | -68447456.4668981573 | 1.55165895740e+12 | 0.0000448141171460574239 | full data |
| V | A | 0.592766703219904303 | -2.93089043260586597 | nan | nan | 211786.563238 | 0.0000729771273296119009 | full data |
| V | B | 0.592746545538598116 | 66.2275800717690729 | -17710.7763060976454 | nan | 478256245.242 | 0.0000628830069835224090 | full data |
| V | C | 0.592770923234743795 | -381.483870930549311 | 363823.590954631032 | -68447456.4668981573 | 1.55165895740e+12 | 0.0000448141171460574239 | full data |
| either | A | 0.592766703219904303 | -2.93089043260586597 | nan | nan | 211786.563238 | 0.0000729771273296119009 | full data |
| either | B | 0.592746545538598116 | 66.2275800717690729 | -17710.7763060976454 | nan | 478256245.242 | 0.0000628830069835224090 | full data |
| either | C | 0.592770923234743795 | -381.483870930549311 | 363823.590954631032 | -68447456.4668981573 | 1.55165895740e+12 | 0.0000448141171460574239 | full data |
| both | A | 0.592766703219904303 | -2.93089043260586597 | nan | nan | 211786.563238 | 0.0000729771273296119009 | full data |
| both | B | 0.592746545538598116 | 66.2275800717690729 | -17710.7763060976454 | nan | 478256245.242 | 0.0000628830069835224090 | full data |
| both | C | 0.592770923234743795 | -381.483870930549311 | 363823.590954631032 | -68447456.4668981573 | 1.55165895740e+12 | 0.0000448141171460574239 | full data |
| M | A | 0.592446366912214791 | 64.2489554442142050 | nan | nan | 211786.563238 | 0.00307353959494817795 | full data |
| M | B | 0.592375155670618899 | 308.565773270761302 | -62567.0359217943985 | nan | 478256245.242 | 0.00307075408625191901 | full data |
| M | C | 0.592080833551292538 | 5713.97333591436001 | -4668990.99399228740 | 826394763.863151439 | 1.55165895740e+12 | 0.00302421427625550776 | full data |

Model spread of those intercepts (max − min over A/B/C at L_min=16, train_max=256):

| observable | n_models | min pc | max pc | spread (max-min) |
|---|---:|---|---|---|
| H | 3 | 0.592746545538598116 | 0.592770923234743795 | 0.0000243776961456790000 |
| V | 3 | 0.592746545538598116 | 0.592770923234743795 | 0.0000243776961456790000 |
| either | 3 | 0.592746545538598116 | 0.592770923234743795 | 0.0000243776961456790000 |
| both | 3 | 0.592746545538598116 | 0.592770923234743795 | 0.0000243776961456790000 |
| M | 3 | 0.592080833551292538 | 0.592446366912214791 | 0.000365533360922253000 |

That spread is **not** a statistical confidence interval.

## Amplitude sign test (Model B/C)

Signs of fitted `a4` and `a6` on every successful training window.
Cancellation is only meaningful if the signed amplitudes are stable
across observable / window / topology rather than accidental averages.

| observable | mode | Model B | Model C |
|---|---|---|---|
| H | shared | n=28; a4 signs {'+': 11, '-': 17}; a6 signs {'+': 17, '-': 11}; median a4/a6=-0.00133896 | n=21; a4 signs {'+': 4, '-': 17}; a6 signs {'+': 14, '-': 7}; median a4/a6=-0.000749446 |
| H | independent | n=3; a4 signs {'+': 2, '-': 1}; a6 signs {'+': 1, '-': 2}; median a4/a6=-0.000943433 | n=1; a4 signs {'-': 1}; a6 signs {'+': 1}; median a4/a6=-0.000193343 |
| V | shared | n=28; a4 signs {'+': 11, '-': 17}; a6 signs {'+': 17, '-': 11}; median a4/a6=-0.00133896 | n=21; a4 signs {'+': 4, '-': 17}; a6 signs {'+': 14, '-': 7}; median a4/a6=-0.000749446 |
| V | independent | n=3; a4 signs {'+': 1, '-': 2}; a6 signs {'+': 2, '-': 1}; median a4/a6=-0.000885967 | n=1; a4 signs {'+': 1}; a6 signs {'-': 1}; median a4/a6=-0.000197630 |
| either | shared | n=28; a4 signs {'+': 11, '-': 17}; a6 signs {'+': 17, '-': 11}; median a4/a6=-0.00133896 | n=21; a4 signs {'+': 4, '-': 17}; a6 signs {'+': 14, '-': 7}; median a4/a6=-0.000749446 |
| either | independent | n=3; a4 signs {'+': 1, '-': 2}; a6 signs {'+': 2, '-': 1}; median a4/a6=-0.000988155 | n=1; a4 signs {'-': 1}; a6 signs {'+': 1}; median a4/a6=-0.000201402 |
| both | shared | n=28; a4 signs {'+': 11, '-': 17}; a6 signs {'+': 17, '-': 11}; median a4/a6=-0.00133896 | n=21; a4 signs {'+': 4, '-': 17}; a6 signs {'+': 14, '-': 7}; median a4/a6=-0.000749446 |
| both | independent | n=3; a4 signs {'+': 3}; a6 signs {'-': 3}; median a4/a6=-0.000963657 | n=1; a4 signs {'+': 1}; a6 signs {'-': 1}; median a4/a6=-0.000209882 |
| M | shared | n=28; a4 signs {'+': 16, '-': 12}; a6 signs {'+': 12, '-': 16}; median a4/a6=-0.00102520 | n=21; a4 signs {'+': 15, '-': 6}; a6 signs {'+': 6, '-': 15}; median a4/a6=-0.000509348 |
| M | independent | n=3; a4 signs {'+': 2, '-': 1}; a6 signs {'+': 1, '-': 2}; median a4/a6=-0.000921807 | n=1; a4 signs {'+': 1}; a6 signs {'-': 1}; median a4/a6=-0.000203217 |

## Noise amplification

L1 = sum |w_i| (fully correlated bound). L2 = sqrt(sum w_i^2) (independent
Gaussian bound). Empirical = batch SD of the weighted combination divided
by batch SD of the largest participating size. Because different L are
independent campaigns, matching covariance lives inside each p_L and the
empirical factor should track L2 more closely than L1.

| method | observable | median L1 | median L2 | median empirical | max empirical | median cond |
|---|---|---:|---:|---:|---:|---:|
| L4_cancel | H | 1.70901 | 1.40228 | 1.27178 | 1.66674 | 2.86724e+7 |
| L4_cancel | V | 1.70901 | 1.40228 | 1.27178 | 1.66674 | 2.86724e+7 |
| L4_cancel | either | 1.70901 | 1.40228 | 1.27178 | 1.66674 | 2.86724e+7 |
| L4_cancel | both | 1.70901 | 1.40228 | 1.27178 | 1.66674 | 2.86724e+7 |
| L4_cancel | M | 1.70901 | 1.40228 | 1.38522 | 1.64377 | 2.86724e+7 |
| L4_L6_cancel | H | 2.43645 | 1.82879 | 1.46739 | 2.01973 | 2.57444e+11 |
| L4_L6_cancel | V | 2.43645 | 1.82879 | 1.46739 | 2.01973 | 2.57444e+11 |
| L4_L6_cancel | either | 2.43645 | 1.82879 | 1.46739 | 2.01973 | 2.57444e+11 |
| L4_L6_cancel | both | 2.43645 | 1.82879 | 1.46739 | 2.01973 | 2.57444e+11 |
| L4_L6_cancel | M | 2.43645 | 1.82879 | 1.63144 | 1.96927 | 2.57444e+11 |
| annihilator_13/4 | H | 2.01198 | n/a | n/a | n/a | 1.49179e+6 |
| annihilator_13/4 | V | 2.01198 | n/a | n/a | n/a | 1.49179e+6 |
| annihilator_13/4 | either | 2.01198 | n/a | n/a | n/a | 1.49179e+6 |
| annihilator_13/4 | both | 2.01198 | n/a | n/a | n/a | 1.49179e+6 |
| annihilator_13/4 | M | 2.01198 | n/a | n/a | n/a | 1.49179e+6 |

## Notes

- Shared wrapping H/V/either/both roots are identical at printing precision,
  so their ordinary/cancellation tables repeat. Cluster matching M does not.
- Independent-mode sizes are 32, 64, 128, 256 only; they appear in the CSVs
  but not in the specified rolling folds.
- Issue #7 L=2..5 totals were read and are the calibration source of Issue #9;
  they are not used as rolling-fit points.
- Fitted `a4` and `a6` signs flip across training windows for wrapping
  observables. The implied correction `a4 L^{-4}` is smaller than the
  batch SE of `p_L` (~1e-4), so the signed amplitudes are noise-dominated
  rather than a stable cancellation pattern. Cluster matching M is noisier
  still; several annihilator roots left `[0.590, 0.595]`.
- Batch combinations skip a batch if any participating size has a non-finite
  root. Empirical noise amp is undefined when fewer than two complete batches
  remain (reported as nan, not 0).
- No hardware or resource recommendation is made.

