# P32 frozen radial and logarithmic challenge

The model challenge used only N=65,85,130 for training and held N=145,170 out
until scoring.  All available independent seed rows were retained.  The
power-correction candidate was selected on training data only; `omega=1` won
over `omega=2,3` by the preregistered training chi-square rule.

| model | fitted parameters | held-out chi-square (3 rows) | improvement over zero |
|---|---|---:|---:|
| fixed `13/8` | `A=0.79768` | 1.058 | 35.27x |
| fixed `13/8` + power | `A=1.1708, B=-3.1689, omega=1` | 1.712 | 21.80x |
| fixed `13/8` + log | `A=0.03687, B=0.17752` | 1.726 | 21.62x |
| fixed `13/8`, H4+H8 | `A4=0.80773, A8=-0.03454` | 1.100 | 33.94x |
| free exponent | `A=0.31200, alpha=1.40615` | 1.661 | 22.47x |

The zero-effect held-out chi-square is 37.32.  The simplest preregistered
`alpha=13/8` H4 law gives the best held-out score and the best conditioning.
The exact rational H8 design column is identifiable and changes sign across
the frozen sizes, but its fitted coefficient is `A8=-0.0345 +/- 0.0542`; adding
it slightly worsens held-out chi-square from 1.058 to 1.100.  The data therefore
support the structural H4 `13/8` model over a zero effect; they do not require
H8, a radial correction, a logarithm, or a free exponent at present precision.

`challenge.json` includes fit covariance, prediction-residual covariance,
training scores, and the held-out zero benchmark.  The three CSV files retain
signed held-out errors, per-size amplitude drift, and nested training-only
leave-one-size-out predictions.
