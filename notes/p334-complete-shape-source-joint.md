# P334: lifetime, full even topology, and the joint-rank carriers

The old-path lifetime signal has a direct full-observer endpoint:
`integral E_top = 1-E W/(N+1)`. It is not confined to a W² contribution inside
J1. The completed same-batch analysis now separates the lifetime mechanism,
marginal thermal shape, connected birth coupling, and the nine joint-rank
carriers without assigning them independent evidence.

## 1. Full E reads the lifetime difference

At N425 the H4-normalized integral E contrast is
`-0.001014130 +/- 0.000344771` on baseline paths and
`-0.000833632 +/- 0.000348650` with the original safe conditional policy.
Their shared uncertainty is retained. These are paired estimators on the same
archive, not a change of population or two independent replications.

The mean-lifetime contrast `0.432019 +/- 0.146872` can be unpacked as

```
H4 Delta K1 = H4 Delta C - H4 Delta W/2 = -0.378430 +/- 0.206010
H4 Delta K2 = H4 Delta C + H4 Delta W/2 =  0.053589 +/- 0.196403.
```

The point direction is predominantly a first-birth displacement, not a large
second-birth delay. Neither endpoint mean contrast separately resolves its
sign; the more precise difference uses their paired covariance. The
squared-lifetime split remains 15.9718 +/- 5.4470 from squared mean and
10.5342 +/- 4.3943 from lifetime variance, as detailed in
`notes/p334-lifetime-square-mechanism.md`.

The conditional policy leaves F1 untouched. Therefore, on every path, safe
minus baseline changes A and E by exactly the same F2 increment. It acts along
the `(1,1)` direction; `E-A=2-2F1` is unchanged. This identity holds at both
p_ref and after integration.

## 2. Where the complete E signal is carried

The nine cells classify the **joint** checkpoint rank `(R_first,R_second)`.
They are not the previous four binary R1-flag states. Their shared indicator
multiplies both orientations, so the constants in A/E cancel within every
cell. All nine mean covariances, including cross-cell terms, are retained.

For N425, safe integral E decomposes into these six transpose groups:

| Joint rank group | H4 contribution +/- shared-batch SE |
|---|---:|
| 00 | -0.000072166 +/- 0.000072658 |
| 11 | -0.000151698 +/- 0.000138928 |
| 22 | 0.000011436 +/- 0.000082599 |
| 01+10 | -0.000434177 +/- 0.000194806 |
| 02+20 | -0.000050214 +/- 0.000097256 |
| 12+21 | -0.000136814 +/- 0.000183270 |

The largest point contribution is 01+10, about half the total. Across these
groups the point contributions mostly reinforce rather than cancel. Within
transpose pairs there is substantial cancellation: the 01 and 10 terms are
`-0.00554991` and `+0.00511574`; the 12 and 21 terms are `+0.00519074` and
`-0.00532756`. The net differences and their errors, not the individual large
signed terms, are the pertinent comparison. These group shares are descriptive
on the already seen archive, not separately selected discoveries.

The 01/10 pairs contain R0 and therefore remain baseline under the existing
whole-pair policy. No attempted conditional replacement of their unknown F1
is hidden in this result. N325 has no resolved full E integral contrast
(`-0.000052446 +/- 0.000413207` safe), with its largest transpose contributions
pointing in opposing directions; no new cross-size exponent is fitted.

At N425 the safe full A/E integral covariance is `-7.09683e-8` (correlation
about -0.2483). The score also includes the complete eight-readout covariance
for baseline/safe, canonical/integrated, A/E. Independent error addition would
discard this structure.

## 3. A lifetime component is not the whole marginal shape or a copula

For the normalized rank mixture Y, energy about fixed p_ref is
`R_ref=E[(Y-p_ref)^2]`, while its ensemble variance subtracts the squared
displacement of its pooled mean. At N425:

| H4 shape coordinate | Mean +/- shared-batch SE |
|---|---:|
| fixed-reference center energy | 2.96770e-6 +/- 2.15135e-5 |
| lifetime energy | 3.65146e-5 +/- 1.21773e-5 |
| total rank R_ref | 3.94822e-5 +/- 2.61371e-5 |
| canonical Beta-mixture Q_ref | 3.95554e-5 +/- 2.61598e-5 |
| ensemble Var Y | 3.89365e-5 +/- 2.58576e-5 |
| mixture mean displacement squared | 5.45733e-7 +/- 6.21790e-7 |

The lifetime and center point terms reinforce; the unresolved total comes
from center uncertainty, not evidence of an opposing center signal. Their
covariance is `3.60158e-11`. Canonical binomial smoothing does not visibly
erase or create the point shape contrast at this precision.

R_ref and Q_ref use the two birth marginals. In contrast,
`center-lifetime` is a joint fixed-reference product, and
`Cov(K1,K2)=Var C-Var W/4` is connected. The latter's normalized N425 H4
contrast is `-1.20899e-5 +/- 2.09091e-5`. Thus the present data do not resolve
a birth-coupling rearrangement or identify a copula from full-A curves.

## 4. The marked source flips do not require a connected winner effect

The source provider decomposes first-birth debt on exactly the accepted R1
population. With r its selection probability,

```
marginal debt M = E(I K1) E(I pi_D) / [r (N+1)]
connected debt B = [E(I K1 pi_D)-E(I K1)E(I pi_D)/r]/(N+1).
```

The collective connected term is exactly `-B`. The four direct/collective
integral point sign flips already occur when B is omitted; restoring B moves
them toward zero. The direct connected H4 contrasts are
`-3.60436e-5 +/- 5.79914e-5` at N325 and
`+1.71018e-5 +/- 4.28858e-5` at N425. Neither resolves a sign. This is a
source-conditioned K1/winner covariance, distinct from the full-population
K1/K2 covariance above. The exact decomposition does not license interpreting
the observed source signs as established.

## Shared source and delivery

Inputs are immutable committed artifacts: lifetime decomposition `be31a113`,
fixed-reference shape `da0080ec`, complete nine-cell A/E `bb79fd47`, and
source-connected decomposition `9ed1e508`. All use the original full births
`9c495ab13e65f2bc93dc0849ee3b73f88724c4b1`; conditional replacements reuse
`0d1e586dafbade5e7d1f9bfc598170d0c881e337` under exactly the old gate.

`results/p334-complete-shape-source-joint/score.json` contains the full 81-cell
mean/SE/covariance, individual cell covariance, six transpose groups, global
A/E table, source/shape/lifetime joint LOO matrix and complete common
covariance. Each nonlinear descriptor uses the provider's original delete-one
batch recentering; no averaging of batchwise conditional products occurs.
Every N retains its twenty original paired 1,000-counter batches. Rank is at
most 19; no inverse covariance or new omnibus test is computed.

Reproduce this committed-vector join with
`/Users/lc/python-envs/research-py311/bin/python scripts/p334_complete_shape_source_joint.py`.
This completes the current handoff with zero new MC, DP, path replay or test
suite, and no additional candidate models.
