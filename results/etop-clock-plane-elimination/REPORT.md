# E_top clock-plane production elimination

## Answer

The one-coordinate `E_top = beta_A A_top` ray is eliminated after the
independent P205 prism is added.  A second measured clock coordinate
changes the model ranking without adding new Monte Carlo samples.

| model | production chi2/df (p) | joint chi2/df (p) | profiled P205 p |
|---|---:|---:|---:|
| A only | 15.5147/7 (0.029939) | 56.4497/10 (1.6908e-08) | 6.7503e-09 |
| A + C | 5.07979/6 (0.53362) | 10.3059/9 (0.32629) | 0.15597 |
| A + W | 5.29222/6 (0.50692) | 16.4874/9 (0.057374) | 0.010716 |
| A + C + W | 4.81181/5 (0.43928) | 9.6447/8 (0.29085) | 0.18445 |

The production-only coefficients are:

- `A+C`: `{'P4_A_top': 0.8287127431200354, 'P4_C': 3.334261982589098}`
- `A+W`: `{'P4_A_top': 0.3171927036592288, 'P4_W': -4.78676029323001}`

Adding `C` to the ray improves the joint profile by `Delta chi2=46.1437` on one df.  Adding `W` also rescues the ray at the declared alpha, but its coefficient-uncertainty-profiled P205 score is near the boundary (`p=0.010716`) whereas `A+C` remains comfortable (`p=0.15597`).  This is common-plane compatibility, not a claim that the production point estimates of the coefficients predict P205: the fixed-coefficient stress rows reject both because the production block leaves those coefficients broad.

Adding `W` after `A+C` gains only `Delta chi2=0.661245` on one df (`p=0.41612`), so the current data do not require a third coordinate.

## Scientific card

- **Mechanism space changed:** universal topology ray versus a topology-plus-measured-clock plane.
- **Result:** one measured clock coordinate absorbs the four-lineage ray failure; after profiling coefficient uncertainty, a common plane remains compatible with the independent quotient prism.  The `C` coordinate is the parsimonious current survivor, while `W` alone is not formally excluded at alpha .01.
- **Not proved:** `C` is not thereby identified with a continuum primary or an independently marked local Potts energy field.
- **Observer / sector / source:** `P4(A_top,E_top,C,W)` / threshold-rank and event-clock plane / aligned production batches.
- **Dependency:** P49/P43/P50/P57 train the coefficients; P205 is an independent raw archive and enters only as a profiled external block.
- **Next discriminator:** record an independently marked same-batch local singlet/energy row and ask whether it replaces `C` and pins the plane coefficients, rather than collecting another untyped E_top ray.

## Reproduction

```bash
python3 scripts/etop_clock_plane_elimination.py --json results/etop-clock-plane-elimination/latest.json --markdown results/etop-clock-plane-elimination/REPORT.md
python3 -m unittest discover -s tests -p 'test_etop_clock_plane_elimination.py'
```
