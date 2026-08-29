# P275 Phase 1 provenance

| block | Huawei DevEnv | N | remote directory | elapsed seconds by i/2i/5i/2 |
|---|---|---:|---|---|
| Zy | DevEnvC_ZyTrST / f415a4bcbd9a438b85f5f29e4a507ea4 | 50 | `/workspace/Matching-One-p275-atop-field-identity` | 108.124 / 109.327 / 104.561 |
| XP | DevEnvC_XPk2PZ / f550f3cb1f774374b6842aa648fda796 | 130 | `/workspace/Matching-One-p275-atop-field-identity` | 215.087 / 215.721 / 217.788 |
| HZ | DevEnvC_HZsCM6 / 033945d8bf8b47a7acf475c595169e07 | 170 | `/workspace/Matching-One-p275-atop-field-identity` | 280.268 / 282.659 / 284.102 |

All hosts checked out runner commit `cb83673fb5f221616a47d53f564635c11e7d0680`
from the same complete git bundle. GCC 10.3.1 produced byte-identical binaries
on the three ARM64 hosts. Each host passed the exact N5 self-test before the
authorization commit was created.

The primary CSV/metadata SHA256 pairs were computed on their originating hosts
and again after transfer. Every value matched; `CHECKSUMS.sha256` records the
local result tree. Metadata validation additionally checked schema, matrix,
Smith invariants `(1,N)`, z representative, seed, counter interval, batch grid,
runner commit, binary hash, root-estimator semantics and same-N priority-field
digests before scoring.
