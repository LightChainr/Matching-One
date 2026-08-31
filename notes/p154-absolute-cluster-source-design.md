# P154: absolute paired-cluster source on the archived Phase-E configurations

Status: estimator design, not a result or a new production assignment.
The primary question is whether a specified **new physical source** changes
the finite directional rank response beyond the declared linear
occupation/Euler controls. It is not a test claiming to identify a new field.

## Source, readout, and the exact limitation

For a configuration on the original graph pair, define

\[
S_\Sigma=c_{\mathrm{black,NN}}+c_{\mathrm{white,matching}},\qquad
s_\Sigma=S_\Sigma/N,\qquad q=I_2-I_0,\quad E=q^2=I_2+I_0.
\]

Count only components containing sites of the corresponding colour; inactive
union-find singletons are not clusters. The white matching graph includes
the NN edges and both square diagonals. The finite positive measure

\[
P_\lambda(\eta)=Z_\lambda^{-1}P_p(\eta)e^{\lambda s_\Sigma(\eta)}
\]

has responses `Jq=Cov(q,sΣ)` and `JE=Cov(E,sΣ)`. Coupling instead to the
extensive `SΣ` multiplies both responses by N: report this normalization,
not a fitted scale power. This is a paired **site-cluster** fugacity, not an
automatic identification with the ordinary one-colour bond-FK Q tangent.
It is matching-even under the full transformation exchanging the two graphs,
colours and `p↔1-p`; this does not assert same-NN-grid complement symmetry
or force Jq to vanish at the square-site reference probability.

The relative source at
[`0912aa43`](https://github.com/LightChainr/Matching-One/blob/0912aa43a9cf32ad9eb3f718a94874efa4ebb3a8/notes/p114-relative-cluster-fugacity-formula.md)
obeys `q=cb-cw-χ`, with `χ=K-T_NN+F4`. Its three-state closure constrains
the **difference**, not the absolute sum. Nevertheless, for these two rank
readouts alone every microscopic source has a conditional projection

\[
\mathbb E[s_\Sigma\mid q]=a+bq+cq^2.
\]

Writing `m_r=E[sΣ|rank=r]`, its coefficients are
`a=m1`, `b=(m2-m0)/2`, `c=(m2+m0)/2-m1`. With `A=E[q]` and `e=E[E]`,

\[
J_q=b(e-A^2)+cA(1-e),\qquad
J_E=bA(1-e)+ce(1-e).
\]

Thus a nonzero response or a two-coordinate rank is not evidence for more
than the existing three topological sectors. The source is specified before
its response is inspected; coefficients varying freely with p, N and
geometry would not constitute a predictive mechanism. Do **not** regress
the source on q/E and then test its q/E covariance: that manufactures zero.

## Exact archived stream, not the earlier one-point pilot

Use the archived mixed-pilot stream whose results are retained at
`e526b9bc5b00cd0e3c17048056792a6ea8a57564`. Its metadata declares
`0578105d92d3822cb48f5c421bd23ff339295cc6`, but that string does not resolve
to the original runner commit and is not treated as a verified freeze pin.
The actual original runner is
`05781051b76001f2b18560d7b0914f2481412584`; its integer-period backend and
the replay backend share blob `22058703c12b168e844088277c9b61d64b9c1d2c`.
Preserve the metadata declaration and this correction separately:

| N | first / second (a,b) | seed | half-open replica interval |
|---:|---|---:|---|
| 65 | (8,1) / (7,4) | 202615465 | [15466000000,15466020000) |
| 130 | (11,3) / (9,7) | 2026154130 | [15466200000,15466220000) |

There are 100 aligned batches of 200 replicas per size and
`p=0.59274605079`. The two sizes are independent blocks; all new source,
rank, control and orientation views within a size reuse its original block.
Inputs are `results/p154-phase-e-mixed-plane-pilot/raw/` files
`n{65,130}_mixed.batches.csv` and `.metadata.json`.

Retain the C++ HNF site numbering and `counter_permutation` from
`src/threshold_rank_integer_period_mc.cpp`. The permutation seed key is
`splitmix64(seed ^ splitmix64(replica + 0xd1b54a32d192ed03))`;
the independent fixed-p count key uses `0x8cb92ba72f3d8dd7` and N comparisons
of the upper 53 RNG bits with p. Both orientations share the permutation
and K. Do not substitute another site labelling, another shuffle, or the
earlier one-point-pilot intervals starting 15465000000/15465200000.

The original B pilot does not store cluster counts or `F4`. The completed
[`c0880c2`](https://github.com/LightChainr/Matching-One/commit/c0880c297b40699563e8be537e777ac8cd4084c8)
Q/R/H replay additionally stores K, K(K−1), occupied NN edges and their
I0/I2 products, but not the new cluster cross moments. A deterministic
reobservation must add the new mark; existing batch means cannot supply it.
No old B/Q/R/H score needs to be rerun.

## Auxiliary linear clock/Euler allocation

Compute the local Euler count independently:

\[
T_{NN}=\sum_{\{i,j\}\in E_{NN}}\eta_i\eta_j,\qquad
F_4=\sum_z\eta_z\eta_{z+e_x}\eta_{z+e_y}\eta_{z+e_x+e_y},\qquad
\chi=K-T_{NN}+F_4.
\]

Each quotient representative anchors one square cell. Matching diagonals
are not additional cells in this Euler formula. `cb-cw-q=χ` is an identity
cross-check, not the definition used to manufacture the control.

For each orientation separately use the fixed control span

\[
Z=\left(K/N,\ K(K-1)/[N(N-1)],\ T_{NN}/(2N),\ \chi/N\right),
\quad \beta=\operatorname{Cov}(Z,Z)^+\operatorname{Cov}(Z,s_\Sigma).
\]

The auxiliary source is
`s_res=sΣ-E[sΣ]-βᵀ(Z-E[Z])`. Report raw and residual responses and

\[
J_{O,\mathrm{res}}=\operatorname{Cov}(O,s_\Sigma)
-\operatorname{Cov}(O,Z)\beta,\qquad
V_{\mathrm{res}}=\operatorname{Var}(s_\Sigma)
-\operatorname{Cov}(s_\Sigma,Z)\beta.
\]

This is a retrospective **linear clock/Euler residual**, not exact
conditioning on K, exact thermal orthogonality, or an RG projection. Even
though the controls include the polynomial span used by Q/R/H, arbitrary
functions of K have not been removed. At these injective square geometries
`E[χ|K=k]=k−2N(k)₂/(N)₂+N(k)₄/(N)₄`; subtracting this known Euler mean does
not provide the unknown conditional mean of the absolute cluster count.

## Minimal data and estimation contract

Store per batch/orientation the first moments and symmetric Gram matrix of
`(q,E,K,K(K−1),T_NN,SΣ,χ)`, plus sample count and the old K1/K2 totals for
stream provenance. This seven-coordinate matrix is sufficient for the
stated estimands. The replay retains cb and cw as well, giving a full
nine-coordinate Gram, and stores F4's first moment separately. No F4 Gram
is missing: every required F4 product follows exactly from
`F4=χ−K+T_NN` and that Gram. This supports the explicit black/white/Euler
decomposition without redundant second moments. `I0*SΣ`
and `I2*SΣ` follow from `(E*SΣ−q*SΣ)/2` and `(E*SΣ+q*SΣ)/2`.

Pool sample sums before computing sample covariances, using denominator
`n_samples−1`. Project the response, not the per-configuration graph:
`P4[J]=(J_first−J_second)/Δcos4`, with exact rational
`cos4(a,b)=(a⁴−6a²b²+b⁴)/(a²+b²)²`. Keep raw orientation responses too.
Delete one whole aligned batch in both orientations, recompute pooled
covariances and each orientation's β, then retain the full joint covariance
and all 100 delete-one vectors. Do not treat the raw/residual or q/E rows
as independent evidence. Residual variance is a property of this declared
linear span, not proof of an independent continuum state.

## What already exists, and what this can decide

P34 at `80fbdd1e9a380a87a3c56dec7795ceebb0ada23e` already stores an
18-coordinate Euler Gram including q, cb and cw in
`results/server-20260828/P34/N{65,85}/mc.motifs.jsonl`.
P40 at `291854a518b03eef4293431b89254f0f4429da53` has a 20-coordinate Gram
and cross-geometry Gram in
`results/local-20260830/P40-production-motif-projection/N{65,85}/mc.motifs.jsonl`.
These million-sample archives can directly supply `Cov(q,SΣ)` and source
variance; they do **not** store `q²*SΣ`, and their recorded `E` means the
occupied-edge count, not E_top. Neither is an N65/N130 norm-4 source analysis.

A resolved residual Jq/JE would show that this specified source couples to
the directional rank readout beyond the declared linear occupation/Euler
span. A raw response that disappears after projection would instead
localize it to that span at this design and precision. An unresolved result
at 20k distinguishes neither exact zero nor absence of other singlets.
No numerical power or response has been inferred in this design note.

This two-size fixed-p response cannot settle the original norm-4 identity.
That comparison still needs the same source mixed moments on the relevant
dyadic generations/second lineage and their joint response uncertainties;
the old K1/K2 histograms do not contain them. Reconstructing a full thermal
response additionally needs per-K source and source×rank profiles (or the
corresponding higher mixed thermal moments). For example,
`∂p Cov(O,sΣ)=κ(O,sΣ,K)/[p(1−p)]` requires a third mixed moment, not just
the seven-coordinate Gram above. Those are explicit further readouts, not
prerequisites for reporting the present new-source response.
