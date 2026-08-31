# Fresh F4 → global U: fixed-budget transmission

Zero-projection family: **NOT_EXCLUDED**.
Finite-resolution action: **INCONCLUSIVE_STOP_FIXED_BLOCK_WITHOUT_TOP_UP**.

| N | V_F4 | paired jackknife SE | simultaneous family95 interval | resolution |
|---:|---:|---:|---:|---|
| 65 | 0.0648901187 | 0.243070401 | [-0.542228153, 0.672008391] | inconclusive |
| 85 | 0.808540693 | 0.381515717 | [-0.144373203, 1.76145459] | inconclusive |
| 130 | 0.0471853018 | 1.36825773 | [-3.37031952, 3.46469012] | inconclusive |
| 170 | -0.735272655 | 2.23245824 | [-6.31129582, 4.84075051] | inconclusive |

## Definition and inference

The sole source is the unscaled number F4 of fully occupied elementary faces, with source measure exp(t F4). Ordinary degree is N; forced-face degree is N−4, including births at forced k=0. J=N p^4(forced−ordinary).
Every central estimate and each of the 100 paired batch deletions refits its own fresh ordinary pooled root inside [0.55,0.65], then differentiates the root, numerator and denominator of U. No old anchor or old sample enters.
Four primary V coordinates use Bonferroni normal critical value 2.49770547441. Independent N seed domains form separate covariance blocks; matching batch numbers across N are not pairing.
The ±0.5 practical band uses the same bulk-source/global-U units. A resolved small nonzero effect can reject zero while all intervals remain inside that band. A non-excluding interval does not establish zero; an overlapping interval is inconclusive.

## Raw/source diagnostics

score.json retains the full central vector, all original delete-one vectors and their covariance. Its three V terms are fixed-p source E jet, moving-root E curvature and denominator response; their errors are correlated and are not added independently.

## Scope

Prediction scope is zero projection at the four declared finite sizes. This is not field identification, an asymptotic model, a source winner, or a validation of another source. No sign was predicted. The fixed block ends here regardless of outcome: no sample top-up, source substitution or size substitution is authorized by this score.

Freeze commit: `0f7a083770d31095e7b4d688d544637d8fc09658`; old_data_pooled: false.
