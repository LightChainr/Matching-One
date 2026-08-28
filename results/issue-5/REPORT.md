# Issue #5: blind finite-size extrapolation audit

This report is a computational record. It does not claim a new
threshold value and it does not attach a statistical confidence
interval to any intercept.

Four quantities are kept separate throughout:

1. **Arithmetic precision** — how far 60/100/160 dps mpmath
   linear algebra agrees on the same fit.
2. **Out-of-sample prediction error** — rolling-origin and
   withheld-tail errors on cylinder widths that were not used
   to fit the model.
3. **Intercept drift** — movement of the fitted infinite-size
   intercept when `n_min` or the training window changes.
4. **Model-to-model spread** — disagreement among correction
   bases. Min/max of an exploratory ensemble is labelled
   **model spread / exploratory range**, never a CI.

Fit residuals are not treated as uncertainties on `p_c`.

## Inputs and protocol

- Data: `data/jacobsen_2015_square_site_cylinder.csv` (n = 1..21).
- Estimator: `scripts/finite_size_audit.py` (unmodified math).
- Models compared in every grid job: `4` / `4,6` / `4,6,8` /
  `4,6,8,10` / `4,6,8,10,12`.
- Grid: dps ∈ {60, 100, 160} × min_train ∈ {5..10} × holdout ∈ {2,3,4}
  → 54 jobs, 270 model summaries,
  1755 folds.
- Runner: `scripts/run_issue5_grid.py` with
  `ProcessPoolExecutor(max_workers=8)`.
- Rolling-origin leakage check (`train_max < test_min` on every fold): PASS.
- No model was skipped on the full n=1..21 grid.

## Baseline reproduction

Command (unchanged): `min_train=8`, `holdout=2`, `dps=100`,
models `4 4,6 4,6,8 4,6,8,10 4,6,8,10,12`.

| rank | model | folds | median RMSE | intercept range | full-fit intercept | score |
|---:|---|---:|---|---|---|---|
| 1 | `4,6,8,10,12` | 5 | 9.32034738575640447528858712194e-12 | 9.54952411024743892590166980115e-11 | 0.592746050900176725024718949109 | -21.0505861708766 |
| 2 | `4,6,8,10` | 6 | 2.23187149787359196519354579520e-11 | 1.60513293599243706984168447466e-10 | 0.592746050975478220449915834470 | -20.4458198077358 |
| 3 | `4,6,8` | 7 | 8.13763548277674795788817005325e-10 | 0.00000000171065607880526133909969151804 | 0.592746052559226703672106027284 | -17.8563390628451 |
| 4 | `4,6` | 8 | 0.00000000224671577277305796560974821216 | 0.00000000194505236909764704115985822869 | 0.592746049501568461679083656558 | -17.3595205668986 |
| 5 | `4` | 9 | 0.000000161009249582943267564105193112 | 0.000000139474218540816927968095592102 | 0.592746195641691671030297139559 | -13.6486552373247 |

Ranking uses the script's diagnostic score
`log10(median_rmse) + log10(intercept_range + floor)`.
That score is not a likelihood and not a confidence statement.

## 1. Arithmetic precision

This section asks only whether the same least-squares problem
changes when mpmath working precision is 60, 100, or 160 dps.
It is **not** a statement about the statistical accuracy of `p_c`.

Compared 90 matched (min_train, holdout, model)
triples present at all three dps values.

The 30-significant-digit JSON serialization written by `finite_size_audit.py` is identical across 60/100/160 dps for 90/90 matched triples. That is the precision of the stored audit tables, not of `p_c`.

Recomputing the same three quantities from the unmodified
`fit_linear` / `rolling_folds` primitives, then comparing the
unrounded mpmath values, shows where the linear algebra itself
stops changing:

| quantity | 60 vs 100 dps | 100 vs 160 dps | 60 vs 100 agree ≥40 decimals | 100 vs 160 agree ≥40 decimals |
|---|---|---|---:|---:|
| full_fit_intercept | worst 55 decimals, median 55 decimals, best 56 decimals | matched at working precision (no difference in printed mpmath values) | 90/90 | 90/90 |
| median_rmse | worst 60 decimals, median 61 decimals, best 63 decimals | matched at working precision (no difference in printed mpmath values) | 90/90 | 90/90 |
| intercept_range | worst 60 decimals, median 63 decimals, best 64 decimals | matched at working precision (no difference in printed mpmath values) | 90/90 | 90/90 |

For `full_fit_intercept`, raising working precision from 60 dps to 100 dps last changes a digit around decimal place 55 in the worst matched triple (median 55).
Raising 100 dps to 160 dps matched at working precision (no difference in the printed mpmath values of `full_fit_intercept`).

Any residual 100-vs-160 dps movement is at working precision
and is far below the out-of-sample prediction error.
Arithmetic noise is not the dominant uncertainty on this sequence.

See `results/issue-5/precision_stability.csv`.

## 2. Out-of-sample prediction error

Rolling-origin folds from the grid (dps = 100 slice):

| model | n configs | median of median RMSE | worst median RMSE |
|---|---:|---|---|
| `4` | 18 | 0.000000273527627593021897172796297593 | 0.00000274282591456471168032388229310 |
| `4,6` | 18 | 0.00000000655835178183997106307075264452 | 0.0000000638683369387251543468797572319 |
| `4,6,8` | 18 | 0.00000000158995343897363655806890601732 | 0.0000000468253419433153329448944808822 |
| `4,6,8,10` | 18 | 1.47480246003842481018553255055e-10 | 0.0000000228261900365097063233154523865 |
| `4,6,8,10,12` | 18 | 3.70677364276521756090111692208e-11 | 0.00000000417073215719298262506297832043 |

Lowest typical rolling RMSE (dps=100): `4,6,8,10,12` (median of median RMSE = 3.70677364276521756090111692208e-11).
Highest typical rolling RMSE (dps=100): `4` (median of median RMSE = 0.000000273527627593021897172796297593).

These RMSE values are prediction errors on withheld widths,
not standard errors of the intercept.

### Blind final-tail experiments

Selection used only training widths. The tail was scored after
the configuration was frozen.

| experiment | train n ≤ | test n | selected model | selected n_min | training-only median RMSE | tail RMSE | tail max abs |
|---|---:|---|---|---:|---|---|---|
| H2 | 19 | 20,21 | `4,6,8,10,12` | 10 | 1.93679646973007062745629538359e-12 | 2.35817630380513701253632319314e-12 | 3.05998960144521601153981727686e-12 |
| H3 | 18 | 19,20,21 | `4,6,8,10,12` | 8 | 1.40536544492378601812717082526e-11 | 1.37758612509855313720848457635e-11 | 1.93748706392088330942548115319e-11 |
| H4 | 17 | 18,19,20,21 | `4,6,8,10` | 6 | 1.68194143257641766368221106045e-10 | 3.98611308606144913729548494078e-10 | 5.26641545853734112703402450882e-10 |

#### H2

All exponent sets, n_min values, model orders, and ensemble weights are chosen from training-only rolling-origin prediction error (holdout=2) on n <= 19. Withheld tail values are read only after the selected configuration and ensemble are frozen.

Selected `4,6,8,10,12` with `n_min=10` (training-only diagnostic score -22.0790120259091).

| n | predicted | true withheld | signed error | absolute error |
|---:|---|---|---|---|
| 20 | 0.592744227386246029164134023877 | 0.592744227384919961820933361318 | 1.32606734320066255887066960687e-12 | 1.32606734320066255887066960687e-12 |
| 21 | 0.592744551484431471604180736475 | 0.592744551481371482002735520463 | 3.05998960144521601153981727686e-12 | 3.05998960144521601153981727686e-12 |

Exploratory top-5 ensemble (training-only; equal-weight tail RMSE 8.20889056325791484252461576305e-12; inverse-RMSE tail RMSE 5.27340245162955672361099655954e-12).

| n | equal-weight predicted | signed error | abs error |
|---:|---|---|---|
| 20 | 0.592744227390468018037115290456 | 5.54805621618192913804928566806e-12 | 5.54805621618192913804928566806e-12 |
| 21 | 0.592744551491569071952384150455 | 1.01975899496486299921755804585e-11 | 1.01975899496486299921755804585e-11 |

0 candidate (model, n_min) pairs were skipped for too few training points and were not filled in.

#### H3

All exponent sets, n_min values, model orders, and ensemble weights are chosen from training-only rolling-origin prediction error (holdout=3) on n <= 18. Withheld tail values are read only after the selected configuration and ensemble are frozen.

Selected `4,6,8,10,12` with `n_min=8` (training-only diagnostic score -20.8120773522100).

| n | predicted | true withheld | signed error | absolute error |
|---:|---|---|---|---|
| 19 | 0.592743810737677870184664739070 | 0.592743810731291551793346908544 | 6.38631839131783052577665663207e-12 | 6.38631839131783052577665663207e-12 |
| 20 | 0.592744227397295436952735497048 | 0.592744227384919961820933361318 | 1.23754751318021357292753812707e-11 | 1.23754751318021357292753812707e-11 |
| 21 | 0.592744551500746352641944353557 | 0.592744551481371482002735520463 | 1.93748706392088330942548115319e-11 | 1.93748706392088330942548115319e-11 |

Exploratory top-5 ensemble (training-only; equal-weight tail RMSE 3.34046759220444018430053231754e-11; inverse-RMSE tail RMSE 3.42883949072973489665096274847e-11).

| n | equal-weight predicted | signed error | abs error |
|---:|---|---|---|
| 19 | 0.592743810749526333152276409315 | 1.82347813589295007708736064167e-11 | 1.82347813589295007708736064167e-11 |
| 20 | 0.592744227416139592423233673733 | 3.12196306023003124146591667162e-11 | 3.12196306023003124146591667162e-11 |
| 21 | 0.592744551526542761973734682386 | 4.51712799709991619225528278231e-11 | 4.51712799709991619225528278231e-11 |

3 candidate (model, n_min) pairs were skipped for too few training points and were not filled in.

#### H4

All exponent sets, n_min values, model orders, and ensemble weights are chosen from training-only rolling-origin prediction error (holdout=4) on n <= 17. Withheld tail values are read only after the selected configuration and ensemble are frozen.

Selected `4,6,8,10` with `n_min=6` (training-only diagnostic score -18.7004322031237).

| n | predicted | true withheld | signed error | absolute error |
|---:|---|---|---|---|
| 18 | 0.592743268120550651488277667054 | 0.592743267887634361790309542529 | 2.32916289697968124524899413773e-10 | 2.32916289697968124524899413773e-10 |
| 19 | 0.592743811068735771264330850234 | 0.592743810731291551793346908544 | 3.37444219470983941690070861089e-10 | 3.37444219470983941690070861089e-10 |
| 20 | 0.592744227820917650248926144891 | 0.592744227384919961820933361318 | 4.35997688427992783572604001684e-10 | 4.35997688427992783572604001684e-10 |
| 21 | 0.592744552008013027856469633166 | 0.592744551481371482002735520463 | 5.26641545853734112703402450882e-10 | 5.26641545853734112703402450882e-10 |

Exploratory top-5 ensemble (training-only; equal-weight tail RMSE 3.90508358053680474671745572660e-10; inverse-RMSE tail RMSE 3.65205028098832030199501625585e-10).

| n | equal-weight predicted | signed error | abs error |
|---:|---|---|---|
| 18 | 0.592743268100737539468895436306 | 2.13103177678585893776740543184e-10 | 2.13103177678585893776740543184e-10 |
| 19 | 0.592743811054085489098319490619 | 3.22793937304972582074544990616e-10 | 3.22793937304972582074544990616e-10 |
| 20 | 0.592744227813089508675905501498 | 4.28169546854972140179465191611e-10 | 4.28169546854972140179465191611e-10 |
| 21 | 0.592744552007725982019194338744 | 5.26354500016458818281016713684e-10 | 5.26354500016458818281016713684e-10 |

10 candidate (model, n_min) pairs were skipped for too few training points and were not filled in.

See `final_holdout_h2.json`, `final_holdout_h3.json`,
`final_holdout_h4.json`, and `fold_errors.csv`.

## 3. Intercept drift

Intercept range is the span of full-subset intercepts obtained
by raising `n_min` inside a fixed model, as implemented by
`intercepts_by_nmin` in `finite_size_audit.py`. It measures
sensitivity to the lower cutoff, not sampling error.

| model | n configs (dps=100) | median intercept range | max intercept range | full-fit spread across configs |
|---|---:|---|---|---|
| `4` | 18 | 0.000000202139369862742010672831800081 | 0.00000144106324511479860914143671616 | 0.00000139092964976649358551789600000 |
| `4,6` | 18 | 0.00000000502485373958581892581262856074 | 0.0000000473107068129398360259730746925 | 0.0000000473107068129398360259740000000 |
| `4,6,8` | 18 | 0.00000000251658329486179035669405212384 | 0.0000000293452897502463322967706252911 | 0.0000000293452897502463322967710000000 |
| `4,6,8,10` | 18 | 3.98395995467984573589057464098e-10 | 0.0000000167504819702663374270984689624 | 0.0000000166698446634713616968650000000 |
| `4,6,8,10,12` | 18 | 1.18829314093560344745071837864e-10 | 0.00000000380657573907357480125212277950 | 0.00000000380657573907357480125200000000 |

Smallest typical intercept range: `4,6,8,10,12` (median range 1.18829314093560344745071837864e-10).
Largest typical intercept range: `4` (median range 0.000000202139369862742010672831800081).

See `intercept_stability.csv`.

## 4. Model-to-model spread

Different correction bases fitted to the same cylinder sequence
do not share one intercept. The spread below is an exploratory
range over models, **not** a statistical confidence interval.

### Top 5 at baseline protocol (dps=100, min_train=8, holdout=2)

Members chosen by training-only rolling-origin median RMSE.
Equal weights. No withheld tail entered the selection.

| model | min_train | holdout | training median RMSE | full-fit intercept |
|---|---:|---:|---|---|
| `4,6,8,10,12` | 8 | 2 | 9.32034738575640447528858712194e-12 | 0.592746050900176725024718949109 |
| `4,6,8,10` | 8 | 2 | 2.23187149787359196519354579520e-11 | 0.592746050975478220449915834470 |
| `4,6,8` | 8 | 2 | 8.13763548277674795788817005325e-10 | 0.592746052559226703672106027284 |
| `4,6` | 8 | 2 | 0.00000000224671577277305796560974821216 | 0.592746049501568461679083656558 |
| `4` | 8 | 2 | 0.000000161009249582943267564105193112 | 0.592746195641691671030297139559 |

- ensemble mean: `0.592746079915628356371224321396`
- ensemble median: `0.592746050975478220449915834470`
- minimum intercept: `0.592746049501568461679083656558`
- maximum intercept: `0.592746195641691671030297139559`
- model spread / exploratory range: `0.000000146140123209351213483001000000`

### Top 5 across the full dps=100 grid by training-only median RMSE

Members chosen by training-only rolling-origin median RMSE.
Equal weights. No withheld tail entered the selection.

| model | min_train | holdout | training median RMSE | full-fit intercept |
|---|---:|---:|---|---|
| `4,6,8,10,12` | 10 | 2 | 2.15316152543163063418936713737e-12 | 0.592746050847292954520075612872 |
| `4,6,8,10,12` | 9 | 2 | 3.41223366125780473357950908890e-12 | 0.592746050865197948813053234549 |
| `4,6,8,10,12` | 10 | 3 | 3.43193863102123952104113654307e-12 | 0.592746050847292954520075612872 |
| `4,6,8,10,12` | 9 | 3 | 5.10162848105311365687872364334e-12 | 0.592746050865197948813053234549 |
| `4,6,8,10,12` | 10 | 4 | 5.11098607454831040621666164887e-12 | 0.592746050847292954520075612872 |

- ensemble mean: `0.592746050854454952237266661543`
- ensemble median: `0.592746050847292954520075612872`
- minimum intercept: `0.592746050847292954520075612872`
- maximum intercept: `0.592746050865197948813053234549`
- model spread / exploratory range: `1.79049942929776216770000000000e-11`

Across every dps=100 grid summary:

- minimum full-fit intercept: `0.592746004135914017836894671744`
- maximum full-fit intercept: `0.592747497230718245011978312900`
- model spread / exploratory range: `0.00000149309480422717508364115600000`

## Files

| path | contents |
|---|---|
| `results/issue-5/environment.txt` | host, python, packages, git HEAD |
| `results/issue-5/baseline.json` | unmodified audit, min_train=8, holdout=2, dps=100 |
| `results/issue-5/raw/` | 54 grid JSON payloads |
| `results/issue-5/logs/` | per-job stdout/stderr |
| `results/issue-5/grid_manifest.json` | job success/fail records |
| `results/issue-5/model_grid.csv` | one row per (dps, min_train, holdout, model) |
| `results/issue-5/fold_errors.csv` | every rolling-origin fold |
| `results/issue-5/model_ranking.csv` | diagnostic rank within each job |
| `results/issue-5/intercept_stability.csv` | intercept range and full-fit intercept |
| `results/issue-5/precision_stability.csv` | 60/100/160 dps comparison |
| `results/issue-5/final_holdout_h2.json` | blind n=20,21 |
| `results/issue-5/final_holdout_h3.json` | blind n=19,20,21 |
| `results/issue-5/final_holdout_h4.json` | blind n=18..21 |

## What this does not show

- No ± interval in this report is a statistical CI.
- Training residuals are not intercept uncertainties.
- Ensemble min/max is model spread / exploratory range.
- Arithmetic agreement at 40 decimals is not 40-decimal knowledge
  of the percolation threshold.

