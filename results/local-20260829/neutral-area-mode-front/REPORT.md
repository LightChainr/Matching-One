# Neutral-area Krawtchouk mode-front diagnostic

## Result

The exact covector identity and the exact `N=10` half-filled oracle pass.
The proposed bounded moving-front diagnostic does not: no committed
square-site size reaches either frozen relative-tail tolerance by order 16.

| N | intrinsic p0 | mean rank gap | best R<=16 | best relative tail | relative tail at R=16 | R_0.05 | R_0.10 |
|---:|---:|---:|---:|---:|---:|:---:|:---:|
| 65 | 0.592731127 | 5.548175210 | 0 | 3.43864 | 5.86173e4 | unresolved | unresolved |
| 85 | 0.592761506 | 6.607616795 | 1 | 3.87643 | 5.60491e5 | unresolved | unresolved |
| 130 | 0.592748111 | 8.699175715 | 1 | 4.67262 | 1.51288e7 | unresolved | unresolved |
| 170 | 0.592743126 | 10.338020465 | 0 | 5.24798 | 1.08540e8 | unresolved | unresolved |
| 185 | 0.592747602 | 10.911188908 | 1 | 5.44106 | 2.04019e8 | unresolved | unresolved |
| 265 | 0.592746312 | 13.725339051 | 1 | 6.34076 | 2.57179e9 | unresolved | unresolved |
| 325 | 0.592742671 | 15.630954964 | 0 | 6.90911 | 1.03767e10 | unresolved | unresolved |
| 425 | 0.592747431 | 18.524434659 | 1 | 7.72514 | 7.39444e10 | unresolved | unresolved |

The rapidly increasing tail is not interpreted as a physical high-mode
amplitude.  It shows that the uniform canonical-area functional is badly
conditioned in a Krawtchouk basis localized at the intrinsic center; exact
degree-`N` cancellations are indispensable.

## Exact controls

- Every histogram independently satisfies
  `sum_k q_k/(N+1)=E[Kplus-Kminus]/(N+1)` to the working precision.
- The exhaustive self-matching `N=10` oracle gives area `1/7`, mean gap
  `11/7`, zero odd half-filled modes, and exact reconstruction at order 6.
- A separate arbitrary degree-12 polynomial test reconstructs its integral
  from the complete general-`p0` covector.

## Reproduction

```bash
python3 scripts/analyze_neutral_area_mode_front.py \
  results/server-20260828/P45-root-amplitude/n65.hist.csv \
  results/server-20260828/P45-root-amplitude/n85.hist.csv \
  results/server-20260828/P49-fullcurve-doubling-100m/raw/n130.hist.csv \
  results/server-20260828/P49-fullcurve-doubling-100m/raw/n170.hist.csv \
  results/server-20260828/P43-heldout-fullcurve-500m/raw/n185.hist.csv \
  results/server-20260828/P43-heldout-fullcurve-500m/raw/n265.hist.csv \
  results/server-20260829/P57-norm5-500m/raw/n325_500m.hist.csv \
  results/server-20260829/P57-norm5-500m/raw/n425_500m.hist.csv \
  --max-order 16 --epsilon 0.05 0.10 --dps 80 \
  --output results/local-20260829/neutral-area-mode-front/analysis.json

python3 tests/test_neutral_area_mode_front.py -v
```

The frozen contract is `analysis/neutral_area_mode_front_manifest.yaml`.
