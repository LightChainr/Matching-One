# P275 q2 local-D4 production provenance

- Frozen prereveal commit: `ebfbf40fd5845f53a92c13caaefacbb7d196a575`.
- Runner commit: `6ecb339f93b878389301a1dc978ae4a38c522b5c`.
- ARM64 binary SHA256 on both machines: `12f3fc0daf7709f518f812038d5a991fc532841be01b2e64fdd1c08fe95a5c77`.
- N65 ran on `DevEnvC_ZyTrST` (`f415a4bcbd9a438b85f5f29e4a507ea4`) in 116.313298832 s.
- N130 ran on `DevEnvC_XPk2PZ` (`f550f3cb1f774374b6842aa648fda796`) in 225.540026269 s.
- Both runs used GCC 10.3.1, eight OpenMP threads, one million paired replicas and 100 batches.
- Seeds/counter intervals were `27520260829/[16000000000,16001000000)` and `27520260830/[16010000000,16011000000)`.
- The committed path tables are the microcanonical sufficient statistics used by the scorer. Their local SHA256 values match the remote values in `REMOTE_SPARSE_SHA256.txt`.
- The sparse per-replica marked-birth streams exceed GitHub's single-file limit and remain at the recorded server paths; their byte counts and hashes are retained below.

The production metadata records the exact matrices, Smith invariants, command,
counter domain, compiler, binary hash, and source commit. Both complement audits
are identically zero, and both stderr logs are empty.
