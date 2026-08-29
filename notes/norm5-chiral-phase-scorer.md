# Frozen covariance-aware norm-five chiral scorer

This scorer is frozen before the 200k Huawei response is revealed.  It does
not divide by the noisy minus-hand response for primary inference.  Instead,
for each exact phase \(q_s=((2+i)/(2-i))^s\), it fits one common complex
normalization \(c\):

\[
R_-=c,\qquad R_+=q_s c.
\]

In real coordinates, the four observations are modeled by

\[
\begin{pmatrix}
\Re R_+\\ \Im R_+\\ \Re R_-\\ \Im R_-
\end{pmatrix}
=
\begin{pmatrix}
\Re q_s&-\Im q_s\\
\Im q_s& \Re q_s\\
1&0\\0&1
\end{pmatrix}
\begin{pmatrix}\Re c\\\Im c\end{pmatrix}.
\]

One generalized least-squares fit uses the complete 4x4 covariance.  Four
real observations minus two fitted real parameters gives `chi_square / 2 df`
for H4, H8 and H12.  Relative likelihoods are
`exp[-(chi2-chi2_min)/2]`; their normalized three-model weights are a compact
ranking, not posterior probabilities over every possible mechanism.

The direct complex ratio and wrapped phase residual are retained only as
readable diagnostics.  Re/Im and the two hands are never scored as independent
votes.  The true reflected-pair conjugacy null is reported separately; the
same-parent plus/minus difference is not relabeled as a null.

The scorer contract records SHA-256
`9b381f7ecf651d482bc4cbb2a63d2217893a44feb2007f06302bc45bc32ade9f`
for the frozen production manifest.  A changed run manifest is rejected before
reading the response.

The correlated synthetic H8 oracle returns `chi_square` numerically zero for
H8, ranks H8 first, strongly rejects H4/H12, and reports the exact reflection
null.  Parallel H4 and H12 fixtures are locked in tests.  A separate regression
shows that deleting covariance off-diagonals changes the GLS result.

After Huawei produces the unrevealed result, score it with:

```bash
python3 scripts/score_norm5_chiral_phase.py \
  --input /workspace/Matching-One/results/huawei-20260829/P226-norm5-chiral-fixedp/chiral_response.json \
  --manifest experiments/p226_norm5_chiral_fixedp_production_20260829.json \
  --contract predictions/norm5_chiral_phase_scorer_20260829.json \
  --output /workspace/Matching-One/results/huawei-20260829/P226-norm5-chiral-fixedp/chiral_phase_score.json
```

Rebuild the committed synthetic oracle with:

```bash
python3 scripts/score_norm5_chiral_phase.py \
  --synthetic-target H8 \
  --manifest experiments/p226_norm5_chiral_fixedp_production_20260829.json \
  --contract predictions/norm5_chiral_phase_scorer_20260829.json \
  --output results/synthetic-oracles/norm5_chiral_phase_scorer.json
```

