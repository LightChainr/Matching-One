# Is the closed source only an explicit rank bias?

## The fixed rank-bias-only prediction has an exact answer

Decision: **`fixed_rank_bias_only_U_alias_rejected`** on the fixed N25 Gaussian pair.
The closed source has original-U response +0.126165363414; its
explicit q source contributes +0.0532475351147. Their difference,
+0.0729178282995, is exactly twice the ambient-null graph-cycle
source response. The rational enclosure of that difference/A
excludes zero.
This answers the fixed prediction V_Sstar=V_q, without choosing a new rank
coefficient or altering the source after observing the split.

## A coefficient-fixed cycle source, not a free residual

The configuration identity is
`Sstar=2*beta_null+q-3*K+2*N+2`, with q=ambient rank−1 and
`beta_null=dim ker[H1(occupied NN graph;R) -> H1(torus;R)]`.
The common K source is a Bernoulli-odds reparameterization: its derivative
of root/slope-normalized U is zero. Consequently
`V_Sstar=2*V_beta_null+V_q` in the same bulk exp(tS) units.
The computed V_beta_null is +0.0364589141498.
The alternate graph-cycle basis gives 2*V_beta1=+0.179412898529;
this is one algebraic change of basis, not a second independent experiment.

Ambient-null graph cycles include zero-winding combinations of nontrivial
cycles. Beta_null is not a count of elementary full faces or the cellular
hole count of a filled-cell complex. Its source is a specified finite graph
statistic; no claim of a local continuum operator follows.

## The original observer and root normalization are retained

The two geometries are (5,0) and (4,3), DeltaCos4=1152/625,
`U=A*Y_p/Q_p` at the pooled Q=0 root, A=25^(13/8)/2.
Per geometry the exact rank algebra q^3=q gives
`j_q=E-q_mean^2` and `j_E=q_mean*(1-E)`.
Their p derivatives, source root motion and both denominator corrections
give the same four-term original-U response formula as the parent reader.
All terms and rational enclosures are saved in latest.json.

Only the saved integer coefficient profiles and the committed Sstar response
were consumed. Their q/E/count coefficient arrays agree exactly, establishing
the same finite ensembles and root. The saved root bracket was reused;
there was no enumeration, sampling, production replay or new root search.

## What this separates and what remains open

The claim is one exact finite-observer source separation. These deterministic
N25 calculations are coordinates of the same complete configuration sets,
not independent statistical votes. The axis Z5xZ5 and tilted Z25 quotients
have different Smith classes; no large-N amplitude or continuum mechanism is
identified. The completed independent F4 block's inconclusive decision is
unchanged, as are the P154/P334 source-specific decisions.

This excludes only the explicit unit-coefficient q alias. One scalar response
cannot exclude a post-hoc fitted c*q source. A nonzero beta_null contribution
can still act through its occupancy/rank conditional means; it does not imply
that a source centered within every K/rank becomes visible to global U.

This equilibrium decomposition does not itself compute interior Xi. The later
[execution delivery f5c4a74a](https://github.com/LightChainr/Matching-One/blob/f5c4a74a20bad8589c39e1034cfb209462110dbe/results/p337-endpoint-defect/score/REPORT.md)
now completes that separate calculation: Xi=-10.755718407564073 and
R=U*U_st−U_s*U_t=27.766563581230237 have nonzero rational enclosures.
The fixed source-independent gain model and mixed thermal-only null fail.
That delivered result supersedes the earlier scorer-only snapshot; no duplicate
defect run is the next assignment. This later narrative update does not change
the coefficient reduction or its stored runtime/hash receipt.

## Sources and reproduction

The fixed action comes from 0d19179f6c6c36fdbb34b2d93e35a9d5fe10dad3:
notes/decimation-closed-cluster-gas-action.md. The complete source response is pinned at
ec01768f520e85f1acfd9d3fde9bcf855477254e:results/p337-closed-source-n25/latest.json. The baseline is
b8d043fc493ab6d7f808d0c074571d2fdd8fb60f:results/decimation-plaquette-u.
The exact term values are reported as a short list rather than a trend chart:
there is one finite-pair comparison and no estimated scaling series.

Run `python scripts/analyze_decimation_cycle_rank_split.py --output-dir NEW_DIRECTORY`.
