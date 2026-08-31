# Joint canonical Q activation reaches original U beyond nearest-neighbour contact

**Two fixed finite nulls are excluded.** For the unchanged canonical
`Kreg=K2+K0` and the original N25 direction pair, the complete response
`J2=partial_logQ partial_epsilon² U|Q1,0` is strictly negative. Its
non-nearest-neighbour part is also strictly negative. The next J2 request
in the previous overview is now a delivered result, not preparation.

| Prescribed original-U coefficient | Value |
|---|---:|
| J2, all distinct-site pairs | -0.0055194314248394015 |
| Four NN displacements | -0.001751074454402799 |
| Other nonzero displacements | -0.0037683569704366022 |

All three rational enclosures exclude zero. Values above include the
positive area prefactor `A25=25^(13/8)/2`; the exact rational intervals
for the corresponding quantities divided by A25 are in
[latest.json](../results/regular-pair-joint-u/latest.json).

## The actual mechanism decisions

The first comparison was named in overview commit
`16a6548af2bb090d9cec7b8e9236e31c0199a3f0`. The full response definition
and translation proof were committed at
`7557da5271f85a69ea5426b61ce7e67b94ee8ff2`; the NN/non-NN secondary,
producer and vector scorer were fixed at
`99b58fc18666cfa6d35b96b52bb84c78dec43a55`, before either count file.

1. An unchanged Q1 baseline with first-Q effective log weight purely
   linear/additive in epsilon predicts **J2=0**. The total result excludes
   that global closure for this canonical interaction.
2. A model whose joint activation reaches original U only through the
   four nearest-neighbour displacements predicts **J2_nonNN=0**. The
   second result excludes that more specific finite global closure.

Neither comparison fits a source coefficient, chooses a sign, changes
the marked-site normalization or searches a radius after scoring. The
NN split is geometrically fixed and computed in the same traversal.
It does not exclude all finite-range contact/OPE mechanisms: a non-NN
pair on N25 can still be microscopically close.

## Complete normalization and topology are retained

Each vacant-site summand is `T0+epsilon Kreg/N`; the occupied summand,
original q=r−1 and E=q² are unchanged. With the delivered signed joint
kernel g, the relevant source is

```text
S2(A)=2/N² sum_(unordered distinct vacant x,y) g_xy(A).
```

The [functional derivation](regular-pair-joint-u-functional.md) proves
that the complete J2 is the ordinary one-source U functional acting on
S2. All canonical marked coefficients vanish at Q1; first-Q products
of separate unactivated marks and background-Q corrections disappear.
Source covariance centering, thermal differentiation, induced root
movement and slope normalization remain. The four measured terms are

```text
direct centered  -0.005529496208369856
root movement   -0.0014591233330132293
source slope    +0.000007758539927132672
root slope      +0.0014614295766165514.
```

This is not `Cov(a_x,a_y)` of the one-mark activation derivatives. It
is also not the fixed-distance spatial mean C, nor the conditional
four-line Gram. Their positive values do not predict the sign of the
directional, root-following, normalized q/E response measured here.

For invariant moments, translation reduces S2 to the fixed-origin sum
`G16/(16N)`. Each geometry needs only its `2^24` origin-vacant occupations,
but uses the full `2^25` population normalization. The q/E marginals in
that same traversal obey the exact relation
`N*sum_origin_vacant O=(N-K)*sum_full O` at each K. This is an accounting
identity against the imported baseline, not a separate enumeration.

Adjacent vacant marks share their actual isolated physical edge ID;
this preserves the eight-port component partition. The noncontact
spatial sampler's distinct-singleton convention is not reused for these
adjacent pairs. The [producer note](regular-pair-joint-u-producer.md)
records the original NN/matching component construction and integer units.

## Relation to the new spatial production

Execution's `branch_only`
[a237968f spatial result](https://github.com/LightChainr/Matching-One/blob/a237968f1d7a82d26b46e83c58179dbba7f1a908/notes/regular-pair-spatial-transmission-result.md)
has independently seeded L32/r8 and L64/r16 blocks, 200000 configurations
and 200 batches each. It rejects the finite L64 noncontact zero null
with `C64=6.85546875e-6`; its 99% Monte Carlo interval is
`[5.2033972758e-6,8.5075402242e-6]`. Its 32 pair readouts per configuration
are correlated and averaged before inference.

That experiment establishes spatial continuation of the same specified
colour interaction. The current exact experiment establishes its
joint contribution to original U, including a non-NN part. Different
geometries and observers keep these as two complementary statements,
not two estimates of the same coefficient or a unique continuum field.

## Execution and the next scientific question

The [single run](../results/regular-pair-joint-u/run.json) took about
4.63 seconds including one compilation and the final rational score.
The concurrent exact traversals took 2.77574 seconds for axis and
2.65815 seconds for tilted. The source table was read, with its pinned
SHA256, from `c29d8bce07e6381be0f948b2c0dd1f640e42d395`; no kernel
regeneration was performed. Old q/E counts, original root, slope and
U/A enclosures were imported without an old source rescore or root search.
There were no random samples, cloud operations or scientific tests.

Reproduce in a new output directory:

```bash
/Users/lc/python-envs/research-py311/bin/python \
  scripts/analyze_regular_pair_joint_u.py --output-dir /path/to/new-result
```

These are new prescribed moments of the same finite N25 populations,
not new independent statistical evidence. The next field question is
whether a fixed, explicitly projected, macroscopically separated
response obeys competing size predictions. Another epsilon derivative,
counterterm choice or reinterpretation of these exact counts is not
automatically the next priority. Existing P154/P334/F4 decisions stand.
