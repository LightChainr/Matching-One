# P250 joint T/A common-annihilator rank gate

**Decision:** `joint_rank3_closes_target_only`.

This is a zero-new-sample, exploratory reuse of one pinned branch-only N505 CSV. All T/A, hand, and charge coordinates share the same 160 aligned batches and the single dependency group `p250-projective-leg-N505-80k`; they are not independent evidence.

## Rank gate

| shared recurrence | d1..d4 training | strict d5 holdout | closes at alpha=0.01 | roots |
|---|---|---|---|---|
| rank 2 | chi2=74.7876, df=28, p=3.8841e-06, cov-rank=32 | chi2=26.5815, df=16, p=0.0463756, cov-rank=16 | no | descriptive_only |
| rank 3 | chi2=13.3136, df=10, p=0.206664, cov-rank=16 | chi2=13.2053, df=16, p=0.657693, cov-rank=16 | yes | descriptive_only |

The recurrence is `s(d)=a1 s(d-1)+...+aK s(d-K)`, with one complex coefficient vector shared by all eight sequences. The fit uses only targets within d1..d4. The d5 residual is computed after the coefficients are frozen.

## Coefficients and descriptive roots

- Rank 2: a1=0.961898-0.125094i, a2=-0.188997+0.0703001i.
  Roots (descriptive_only): 0.678586-0.0369036i (complex jackknife SE 0.0793916), 0.283312-0.0881905i (complex jackknife SE 0.137756).
- Rank 3: a1=1.44013+0.373742i, a2=-0.554216-0.154015i, a3=0.0537745-0.0549735i.
  Roots (descriptive_only): 0.792616+0.412544i (complex jackknife SE 1.05923), 0.628083+0.0952714i (complex jackknife SE 0.380896), 0.0194322-0.134073i (complex jackknife SE 0.296331).

## Statistical contract

Each delete-one replicate removes the same batch from every T/A x hand x charge coordinate, refits the common coefficients, and recomputes both training and d5 residuals. The JSON stores the complete jackknife covariance matrices. The covariance-aware quadratic forms use a fixed relative eigenvalue cutoff of 1e-10.

The exact coordinate change is `X=T+A`, `Y=T-A`; therefore common recurrence closure is basis-equivalent in T/A and X/Y. T/A is used only as the stored numerical basis.

## Provenance and boundaries

- Manifest: `analysis/p250_joint_ta_hankel_manifest.json`.
- CSV: `b62b4efd4997bd5f9923949b4e07d8011823c9ed:results/huawei-20260830/P250-z5-projective-leg-cross-scale-n505-80k/response_80k.batches.csv`.
- CSV SHA256: `8b0e06f3fdc577c362e6f2404db60933d0cf489ee532c053288a03c44dc7fe5c`.
- Samples/batches: 80000 / 160 aligned batches.
- Analysis worktree HEAD: `9ae950869becf0e3faf9794f06e20af72314d267`.
- The input remains branch-only; this analysis does not merge or copy its CSV into the Draft branch.
- Rank 4 is intentionally not attempted because five distances would leave no honest holdout.
- A rank-3 pass is finite-sequence closure, not a three-field, C4, or Z5 fusion theorem.
- A rank-3 failure excludes the stated shared scalar rank-at-most-3 recurrence on these   five distances, but does not establish noncommuting directional transfer operators.
- Root labels are descriptive; an ambiguous delete-one assignment is emitted as   `root_not_interpretable` and is never used in the rank decision. Even a matchable   label is not a precise root estimate when its reported jackknife SE is large.
