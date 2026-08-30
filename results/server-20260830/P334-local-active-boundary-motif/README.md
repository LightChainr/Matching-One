# P334 local active-boundary motif pilot

## Frozen question

After the global boundary aggregates in `082dcab` failed, can radius-one
structure localized at the realized first/second homology-completion sites absorb
the age residual left by the fixed cross-size rank-one morphology state?

The preregistered readout has four projective-line-safe counts at each site:
occupied neighbours, essential contacts, same-side essential-contact pairs, and
opposite-side essential-contact pairs.  Their site sum/difference gives four
complement-even and four complement-odd coordinates.  The future exit selects a
site, but its field is evaluated on the frozen current-`k0` state.

## Production

| Size | Host | Samples / batches | At-risk rows | Runtime | Raw SHA-256 |
|---|---|---:|---:|---:|---|
| N325 | HZsCM6 | 20,000 / 20 | 18,290 | 0.276 s | `d1339ec93b9847748a26172beeb6f1fafcbc15ce4ba4fb0cd766bbf1633d90cf` |
| N425 | TgFr7R | 20,000 / 20 | 17,960 | 0.362 s | `cbfec64a52b04b60b8022b63b22428da4cab2def49a3f888a6a812a6015fcbb2` |

Both remote self-tests passed and both stderr logs are empty.  Seeds, counter
ranges, matrices and `k0` are frozen in
`analysis/p334_local_active_boundary_motif_freeze.json`.

## Result

The local motif **fails the common absorption gate**.  Residual retention in
`(N325 first, N325 second, N425 first, N425 second)` order is

`(1.556, 0.823, 0.538, 0.229)`.

The size-joint residual tests give `p=0.0274` at N325 and `p=0.924` at N425.
Only one of four environments falls below the frozen 25% retention target; in
N325-first the motif moves the residual in the wrong direction.  Incremental
within-line R-squared is small but nonzero, `0.0139--0.0177`, so these local
coordinates contain structure without supplying a common cross-size age state.

The frozen 100k extension condition is not met (not all four retentions are below
0.50 and the common gate fails).  No extension was launched.

## Scientific boundary

This closes the minimal radius-one, primitive-line-relative completion-site
candidate, not every possible local history.  Together with `082dcab`, it says
that neither global frontier organization nor this smallest localized
completion-site morphology is the missing common age coordinate.  It does not
prove intrinsic temporal memory and it is not combined with other dependency
groups as independent evidence.
