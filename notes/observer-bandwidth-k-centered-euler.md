# K-centered Euler observer on the square-cell torus

For a periodic square-cell geometry with `N` sites, `2N` distinct nearest-
neighbor edges, and `N` four-distinct-vertex faces, the global Euler observer
is

`O_E = K - E_NN + F_0`.

Uniformly on the `K=k` slice, each edge is fully occupied with probability
`(k)_2/(N)_2` and each face with probability `(k)_4/(N)_4`. Hence

`E[O_E|K=k] = k - 2N (k)_2/(N)_2 + N (k)_4/(N)_4`.

The exact `L=3` oracle enumerates all 512 configurations and verifies this
formula, with zero centered sum on every one of the ten K-slices. It then
checks `p=1/3, 2/5, 1/2`. Translation invariance makes the original degree-one
site projections identical; subtracting the conditional mean removes all nine
of them exactly. Against a deterministic full-cube source, only Walsh degrees
2, 3, and 4 remain.

This is a three-mode kinematic envelope under the declared product-noise clock,
not evidence for three physical fields. Degenerate quotients with repeated face
vertices require their own incidence polynomial and are outside this control.
