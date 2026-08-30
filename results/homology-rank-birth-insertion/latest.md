# Canonical homology-rank-birth insertion

The two elementary gates are typed before any CFT interpretation:

```text
I_01 = 1[r0=0 and r1>=1]
I_12 = 1[r0<=1 and r1=2]
Delta_v r = I_01 + I_12.
```

Thus Russo differentiation gives exactly `M'(p)=f_01(p)+f_12(p)`. A direct
`0->2` jump contributes once to each gate and needs no artificial intermediate state.

| geometry | weighted insertions | transitions | M'(1/2) | f01(1/2) | f12(1/2) |
|---|---:|---|---:|---:|---:|
| axis-L2-degenerate | 32 | 0->0:8, 0->1:8, 0->2:4, 1->2:8, 2->2:4 | 3 | 3/2 | 3/2 |
| gaussian-2-1 | 80 | 0->0:25, 0->1:30, 1->2:20, 2->2:5 | 25/8 | 15/8 | 5/4 |
| axis-L4-fixed-root | 524288 | 0->0:250096, 0->1:80128, 0->2:4624, 1->1:96736, 1->2:45312, 2->2:47392 | 4209/1024 | 5297/2048 | 3121/2048 |

For `0->1`, `ell` is the new primitive rank-one image and `iota` its integral
saturation index. For `1->2`, the canonical mark is the rank-one plateau line
immediately before the second birth. The same `ell` labels both endpoints of every
nonempty essential-H1 interval. A simultaneous `0->2` jump has no canonical
intermediate line and is recorded with `ell=null`.

The axis-L4 control also attaches the existing radius-one landing-sector H4 mark
and the exact spin-4 harmonic of the physical winding line. These remain typed
coordinates of the same insertion, not extra evidence rows.
