# Exact contact structure and source strength, prefix by prefix

This extract supplies named contact descriptors and the **complete-census**
physical score Gram, aligned exactly to the row order of the quartet response
archive `375cd3a12b2b7a87d79148a59f62b95898f9e471`. Its purpose is to let the
subsequent fixed-rankcell analysis distinguish source strength from contact
structure. No response is recalculated and no regression or model decision is
made here.

## Inputs and coordinates

- Identity only (`counter,batch,rankcell,old_rank,k0,d`) from the two committed
  `results/p334-exact-score-quartet-moments/N{N}.npz` files at `375cd3a1`.
- Complete joint-safe vacant-label counts from
  `ac5761ce504c3cd170fa42c86c17d6fb87f0375b`:
  `experiments/p334-finite-source-20260831/census/N{N}/census.csv.gz`.
- Reader: `475738c6`, `scripts/p334_exact_prefix_structure.py`.

Each next-label class `a=(e_first,e_second)` has count `n_a`; `d=N-k0` includes
**all** vacant labels, including labels outside joint safety. On joint-safe
labels let `L_i=1{oldrank_i=0}(e_i-c_i)`. The six structure descriptors are

\[
\sum_a n_a/d,\quad \sum_a(n_a/d)^2,\quad
\sum_a n_a e_{a,f}/d,\quad\sum_a n_a e_{a,s}/d,\quad
\sum_a S_{a,f}/d,\quad\sum_a S_{a,s}/d,
\]

where `S_ai=sum_{u in a} L_i(u)`. The degree and loop entries are **not**
divided again by joint-safe mass. Class collision includes the safe degree
classes only; it does not add the outside-safe category as another class.

## Physical source Gram from exact integer counts

Write `Q_aij=sum_{u in a} L_i(u)L_j(u)` and
`s_Li(u)=(n_a/d)(L_i(u)-S_ai/n_a)`, zero outside joint safety. Then

\[
E_U[s_{Li}s_{Lj}]
=\frac1{d^3}\sum_a\left(n_a^2Q_{aij}-n_aS_{ai}S_{aj}\right).
\]

The implementation saves every numerator as an integer, together with its
`d^3` denominator. It neither estimates these moments from eight quartets nor
forms a Gram of response estimates. Cross entries retain their signs.

For `s_plus=(s_Lf+s_Ls)/2`, `s_minus=(s_Lf-s_Ls)/2`, the stored derived entries
are `(Gff+2Gfs+Gss)/4`, `(Gff-2Gfs+Gss)/4`, and `(Gff-Gss)/4`.
These three display coordinates are exact transforms of the physical Gram,
**not three additional independent predictors**. In cells with only one R0
orientation, only that physical source can be nonzero; the equal plus/minus
energies there follow from this definition, not a fitted symmetry.

## Files and schema

`results/p334-exact-prefix-structure/N325.npz` and `N425.npz` contain:

- `features[20000,12]`, float64, with exact order in `metadata.json`:
  `joint_safe_mass,class_collision,safe_degree_first,safe_degree_second,`
  `safe_loop_first,safe_loop_second,score_gram_ff,score_gram_fs,score_gram_ss,`
  `score_energy_plus,score_energy_minus,score_gram_plus_minus`.
- Original identity arrays in exactly the same row order as `375cd3a1`;
  all nine rank cells and original twenty batches remain intact.
- `score_gram_physical[20000,2,2]` and
  `score_gram_plus_minus[20000,2,2]`.
- Per-class `class_count,class_loop_sum,class_loop_product_sum,class_mass,`
  `class_loop_mean,class_score_gram_numerator`. The class numerator gives its
  contribution to the **unconditional** Gram with denominator `d^3`.
- Aggregated Gram numerators/denominators and integer numerators for the six
  structure descriptors. Empty class means are set to zero and have zero mass.

The two compressed files are 2.973 MB and 3.253 MB. Their generation took 1.38
seconds in the managed Python research environment. Metadata records exact
source/reader commits, input/output SHA256, formulas, axes and dtypes. No fork
gzip was read; no tail, response, finite policy, local determinant, DP, or
fitting task was executed.

Descriptors are exact conditional on each recorded prefix. That does not make
the sampled prefix population exact or create an independent evidence block;
the downstream analysis must keep the shared original batch dependence.
