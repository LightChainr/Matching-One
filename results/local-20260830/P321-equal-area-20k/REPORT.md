# P321 N144 equal-area covariance smoke

- machine: Apple M4, Darwin arm64;
- compiler: Apple clang 21.0.0, `-O3 -DNDEBUG -std=c++17`;
- engine source commit: `65b383084ce0068e722315f5d979320a3fdb221e`;
- sampling: 20,000 replicas per rectangle pair, 20 batches, one thread;
- seed/counter: `32114420260830`, `[0,20000)`;
- shapes: square paired with `rho=16/9,9/4,4,9`;
- engine time: 1.162 seconds for all four pairs;
- covariance-score time: 12.74 seconds;
- repeated-square histogram and moment gates: passed byte-for-byte;
- full root covariance: retained in `score.json`;
- scale/E4 score: not run; N144 alone is insufficient by construction.

Contrast SEs at 20k are `4.76e-4`, `5.91e-4`, `5.65e-4`, and `8.93e-4`.
Contrast correlations range from `0.306` to `0.808`, confirming that the
eventual aspect-ratio curve must use the full aligned covariance.

This is a variance/runtime smoke only.  Root differences were not used to
alter the frozen `N^-2/N^-3` model or conditional E4 prediction.

