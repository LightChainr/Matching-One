# Exact complement/Alexander sign control

This control freezes the discrete sign convention required by Issue 439. A
birth row transforms as

```text
(tau1,tau2) -> (1-tau2,1-tau1),
```

while plateau lines follow a declared bijective involution. At paired,
tie-free thresholds `p` and `1-p`, the exact sectors obey

```text
(P0,P1,P2)_dual = (P2,P1,P0),
M_dual(1-p) = -M(p).
```

With the ordered orientation names and covectors held fixed, the synthetic
H4 contrast changes from `-1/6` to `+1/6`. Closing each orientation under the
same complement transform gives exact midpoint `M=0` in both orientations and
zero odd H4 contrast.

Thresholds equal to a discrete birth time fail closed because inclusive birth
conventions do not by themselves define a complementary tie split. The parser
also rejects incomplete or non-involutive line maps, orientation descriptor
drift, floats, and malformed birth rows.

## Boundary

This is an exact synthetic convention check. It does not validate signs in a
production archive, estimate an H4 amplitude, fit a common radial state or
transfer law, or support a continuum/physics claim. Issue 439 remains open.
