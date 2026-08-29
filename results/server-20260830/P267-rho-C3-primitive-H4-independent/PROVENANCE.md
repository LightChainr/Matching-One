# Provenance

- Branch: `research/p267-rho-c3-primitive-h4-20260829`
- Frozen runner commit: `ac2fb6cd8469793a816d6fc7342034386acd899a`
- Authorization/deployment commit: `a638321e80e5454b1d45053d2cfd50301944e503`
- Transferred git-bundle SHA-256: `f0ce419f968817ef743e9ee581cf0c1cbc880bd1835e124a20940560e1bff668`
- Manifest: `experiments/p267_rho_child_primitive_h4_independent_2m_20260830.json`
- Runtime: Python 3.9.9 on three Huawei Kunpeng DevEnv containers
- Start: `2026-08-29T16:07:38Z`; completion observed: `2026-08-29T16:10:01Z`

| child | machine | PID | seed | counter interval | run SHA-256 | batch SHA-256 |
|---|---|---:|---:|---:|---|---|
| `2omega` | `DevEnvC_ZyTrST` (`f415a4bcbd9a438b85f5f29e4a507ea4`) | 26907 | 2671562001 | [17000000000, 17002000000) | `dd5a5e35a2196a1c21bbcfa76c580fe9b118c29a7800eeea73bf7ee3971c4c8e` | `5de72b1b6f4e5001d5db4ba35009ef06bc8a0844945452db93fe6bde4ddff496` |
| `omega_over_2` | `DevEnvC_XPk2PZ` (`f550f3cb1f774374b6842aa648fda796`) | 16731 | 2671562002 | [17002000000, 17004000000) | `e499aa170353895eaa06667e2ec91ce40e13758ffc20b4163baca5d8339495fe` | `f895c65a7172156fbf004ffefd36d72e380cca858961fe4fa462004ec2b85039` |
| `omega_plus_1_over_2` | `DevEnvC_HZsCM6` (`033945d8bf8b47a7acf475c595169e07`) | 14124 | 2671562003 | [17004000000, 17006000000) | `8b51b30b817ef38cdab454025b422da82308f01b5686684be839a9ec88a57f20` | `4706054733404b9460242190df8539134bd68b811033bfb7bab50c165cd90505` |

Each job used 2,000,000 samples, 100 batches, and 8 workers. Remote and
local SHA-256 values agree. All stderr files are empty and all homology
invariant-failure counts are zero. The distinct seeds and adjacent disjoint
counter intervals make the frozen combined covariance block diagonal across
children; each child retains its measured full real/imaginary 2x2 block.

