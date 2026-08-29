# Matching as a relative FK cluster-fugacity derivative

There is a direct finite partition function for the matching observable. It is not the derivative of one ordinary Potts cluster fugacity; it is a relative black/white cluster source with a fixed Euler counterterm.

For a site configuration `omega` on the square torus, let

- `k_b` be the number of occupied NN clusters;
- `k_w` be the number of vacant matching NN+NNN clusters;
- `L=V-E+F0`, using occupied sites, occupied NN edges, and occupied square faces.

Define the Laurent partition function

```text
Z_rel(Q,p)
 = sum_omega P_p(omega)
   Q^k_b Q^(-k_w) Q^[-(V-E+F0)].
```

The configuration-level Euler/Betti identity gives

```text
q(omega)=k_b-k_w-(V-E+F0),
```

where `q` is also the common primal-minus-matching wrapping charge on the declared finite quotients. Consequently

```text
Z_rel(Q,p)=Q Z_+(p)+Z_0(p)+Q^-1 Z_-(p),
```

with `Z_r` the restricted probability of `q=r`, and exactly

```text
(Q d_Q) log Z_rel |_(Q=1) = E_p[q] = M_N(p),
(Q d_Q)^2 log Z_rel |_(Q=1) = Var_p(q).
```

This is the requested partition/topological-sector derivative formula. It needs no continuum assumption and no fitted coefficient.

On the finite quotients the repository's wrapping differences are configuration-identical. For a continuum modular statement, however, the typed #114 contract still applies: use `cross` or `either`; the finite `direction`/`both` equalities do not make those basis-dependent labels modular scalars.

The local counterterm factorizes explicitly:

```text
Q^[-(V-E+F0)]
 = product_sites Q^-n_x
   product_NN_edges Q^(n_x n_y)
   product_square_faces Q^(-n_1 n_2 n_3 n_4).
```

Thus the source consists of three separately meaningful pieces: black cluster fugacity `+k_b`, complementary white matching-cluster fugacity `-k_w`, and the explicit site/edge/face Euler derivative. At fixed `p`, these are the complete terms. Following a `Q`-dependent critical manifold adds the separate #258 measure-score covariance; continuing a generic-`Q` homology projector may add a projector derivative. Neither should be silently folded into this fixed-`p` identity.

## Why one ordinary FK derivative is insufficient

The committed `N=25`, `(a,b)=(4,3)` Gaussian torus contains two explicit configurations, masks `0x1d24768` and `0xf6aca0`, with

```text
k_black = 2 for both,
the complete histogram of all 16 oriented 2x2 plaquette patterns identical,
q = 0 and -1 respectively,
k_white = 1 and 2 respectively.
```

Therefore no score of the form

```text
a k_black + sum_x f(oriented 2x2 plaquette pattern at x)
```

can equal the matching charge on every configuration, for arbitrary coefficient `a` and arbitrary plaquette-local function `f`. This is a finite exact obstruction to “ordinary one-colour FK derivative plus local plaquette counterterms.” The missing information is precisely the complementary white connectivity or an equivalent topology projector.

This clarifies #233: the relative source is the first `Q`-tangent pull-through obstruction, and complement/matching exchange sends `Q` to `Q^-1`. It also matches #123: composition generates partial boundary partitions, so a paired black/white connectivity state is unavoidable; the all-or-none local tensor cannot remain closed.

For fixed transfer width, tracking both connectivity partitions gives a finite state representation. No claim is made here that its bond dimension remains width-independent, that the continuum relative source is already a known CFT field, or that generic-`Q` projector derivatives vanish.
