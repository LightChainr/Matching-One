# q2 local landing-D4 UV-annihilator

Primary: paired-orientation mean local-D4. `q=A_top` is an exact contact control only.

| size | intrinsic p | P4 D shell | P4 thermal S shell | contact max residual |
|---:|---:|---:|---:|---:|
| 65 | 0.5927073236840990434243671 | 0.06358626822982737049558234 | 0.5623105438662612123432539 | 7.77876909732642713393008e-62 |
| 130 | 0.5928572649994181192332681 | 0.04287211350847562309196797 | 0.1829246731911447759101258 | 1.555753819465285426786016e-61 |

## Heldout q2 score

- target child/parent: `-0.3242098886627524164834385`
- raw P4-shell child/parent: `0.6742354080839887026709753`
- N65-trained thermal beta: `-0.02421340384773446363569725`
- N130 residual: `['0.02254034291492579188279197', '0.1211500483242092653252702']`
- chi-square: `72.86329633810854113924015/2`, survival p `1.506385870771878671531777e-16`
- preregistered decision at alpha `0.01`: **rejected**

## Boundary

- R4-R2 cancels only an R-independent local contact term under the frozen cutoff model.
- The N65 thermal-shell nuisance and H4 amplitude are trained before the N130 child is scored.
- The q-product contact identities are exact controls, not evidence for a field coupling.
- A surviving two-coordinate heldout score nominates the local source; it does not prove Q4 epsilon.
