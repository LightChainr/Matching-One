# All-order topology-insertion no-go

For `q in {-1,0,1}`, every topology-only observer is exactly
`f(q)=a+bq+cq^2`. Insertion defines `S=Delta q`, `D=Delta(q^2)`,
and the gate algebra has basis `(1,S,D,D^2)` with

```text
S D = D,       D^2 = 2S-S^2,
S(S-1)(S-2)=0, D(D-1)(D+1)=0.
```

For any root-occupation-independent mark `H` and insertion polynomial
`G` with `G(0,0)=0`,
write `mu_G=<HG>`, `mu_SG=<HSG>`, `mu_DG=<HDG>`. Then

```text
Cov(f(q),HG) =
 b[(p-1/2)mu_SG + mu_DG/2 - <q>mu_G]
+c[mu_SG/2 + (p-1/2)mu_DG - <q^2>mu_G].
```

This is an all-order no-go: no lookup table, nonlinear function, or
higher polynomial of ambient rank can create an independent insertion
matrix element. The minimal polynomial `Q(Q-1)(Q+1)=0` reduces them all.

The tiny axis/Gaussian oracle performs 8960 exact marked pointwise checks with zero failures.

The escape is precise: use an observer varying within a fixed q sector—
a separated local landing/arm observable, bulk Betti/Euler data, a charged
seam or winding character, or an independent modular/stress response.

## Scientific card

1. MECHANISM SPACE: every q-only connected insertion response collapses to a four-dimensional gate algebra.
2. NOT PROVED: the theorem says nothing about observers that resolve configurations within one q sector.
3. OBSERVER-SECTOR-SOURCE-GEOMETRY: f(q) | ambient-rank topology | arbitrary marked S/D gate polynomial | any finite torus.
4. DEPENDENCY GROUP: previous A_top/J_D and q2 orthogonal scores are corollaries, not independent evidence.
5. UPWEIGHT OBSERVATION: an independent bulk, separated-local, or charged-seam observer with same-field source products.
