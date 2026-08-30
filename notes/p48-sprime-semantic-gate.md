# P48 S-prime scorer semantic gate

## Exact scope

The existing P48 prospective scorer compares `P4[S_prime]` source predictions
with fresh `N=185,265` targets. Both sides are the same typed observable:

```text
channel             cross
combination         even
coordinate          p
orientation order   first_minus_second
normalization       angular_normalized
quantity            orientation_contrast
```

The canonical `map_observable` registry therefore returns exact affine map
`scale=+1`, `offset=0`. The typed entrypoint checks this before importing and
running the frozen numerical kernel.

## Traceability closure

The 2026-08-28 scoring manifest and the separate model artifact are both
declared explicitly by the semantic gate. They and the frozen kernel remain
byte-for-byte unchanged. The canonical audit row now points to the typed
entrypoint rather than implying that the untyped frozen kernel performs the
semantic validation itself.

The typed entrypoint requires `--output`. After the frozen kernel succeeds, it
adds the semantic manifest path and status, both artifact paths, the source and
target descriptors, the applied exact map, and the validation order to that
JSON result. A pre-existing semantic block is rejected rather than overwritten.

## Claim boundary

This is a protocol integration result, not a new numerical or physical result.
It does not score target data, validate the provenance of target samples,
change the frozen model hierarchy, register a new channel identity, or certify
any scorer outside the P48 `S_prime` path. Those repository-wide tasks remain
open under Issue #146.
