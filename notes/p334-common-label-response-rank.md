# Two Euler-invisible controls resolve two mean-response directions

The new shared-label experiment finds a useful distinction. A symmetric
contact tilt has a clear common A response, but no resolved near-reference
H4-normalized difference. An antisymmetric tilt does produce a clear A
difference at both sizes. In the underlying orientation coordinates this
is a nearly diagonal two-control response, not evidence that the symmetric
background was itself a native H4 field.

## Direct common-label responses

Every entry uses one label policy that preserves the full joint immediate
rank/Euler distribution of both geometries. Define g_plus=(L_f+L_s)/2 and
g_minus=(L_f-L_s)/2, with L the R0-only safe-loop mark. Outputs are
S=(Y_f+Y_s)/2 and D=(Y_f-Y_s)/delta_cos4.

| A(p_ref) derivative | N325 +/- original-batch SE | N425 +/- original-batch SE |
|---|---:|---:|
| plus -> S | -5.41618e-5 +/- 7.71650e-6 | -5.87259e-5 +/- 1.26600e-5 |
| plus -> D | +3.67632e-5 +/- 2.69786e-5 | -5.11096e-6 +/- 1.90100e-5 |
| minus -> S | -7.76062e-6 +/- 1.10961e-5 | -3.15439e-6 +/- 1.34342e-5 |
| minus -> D | +1.53657e-4 +/- 3.46721e-5 | +1.48394e-4 +/- 2.42130e-5 |

For E(p_ref), plus -> S is positive at both sizes,
`1.88401e-5 +/- 5.81637e-6` and `1.93281e-5 +/- 5.91970e-6`,
while plus -> D is unresolved. The minus -> D E response is
`-1.93689e-5 +/- 1.84460e-5` at N325 and
`-5.22949e-5 +/- 1.56843e-5` at N425. This is an imposed differential
source response; it does not establish the origin of the unperturbed
global H4 mean.

## Undo the parity coordinates

Let J have rows (future geometry first,second) and columns (input L_f,L_s).
It is the Jacobian of the two mean future responses with respect to the
two shared-label source strengths at zero. Write the parity response
matrix with rows (S,D) and columns (plus,minus) as R. Then

```
J_ff = R_S+ + R_S- + delta*(R_D+ + R_D-)/2,
J_fs = R_S+ - R_S- + delta*(R_D+ - R_D-)/2,
J_sf = R_S+ + R_S- - delta*(R_D+ + R_D-)/2,
J_ss = R_S+ - R_S- - delta*(R_D+ - R_D-)/2,
det J = 2*delta*det R.
```

For A(p_ref), entries of J in units of 1e-4 are

```
N325: [ -1.34611 +/- .27353,  -.01779 +/- .13722 ]
      [  +.10766 +/- .18037,  -.91023 +/- .24368 ]

N425: [ -1.25849 +/- .22196,  +.12961 +/- .20582 ]
      [  +.02088 +/- .17357, -1.24104 +/- .30017 ]
```

Both diagonal effects are negative. The off-diagonal A effects are not
resolved. The full common-label mixture can therefore have a strong
source-odd/observer-difference response while its source-even difference
is weak: each source principally loads its own geometry's future A.
This is a measured pattern, not an exact vanishing-cross-coupling theorem.

## A finite two-control rank witness

The determinants, with delete-one-original-batch propagated errors, are

| Response | N325 det J +/- SE | N425 det J +/- SE |
|---|---:|---:|
| A(p_ref) | 1.22719e-8 +/- 4.07190e-9 | 1.55912e-8 +/- 4.69216e-9 |
| integral A | 4.20205e-10 +/- 1.62825e-10 | 5.14308e-10 +/- 1.28723e-10 |

The point-estimate oriented area of unit-normalized response columns is
.9982/.9927 near p_ref and .9990/.9974 after integration. These angles
are descriptive; symmetric Gaussian intervals near their bounded endpoint
should not be read as exact confidence regions. Input mark correlations
in the same masked covariance are only .0177/.0241, and both input Gram
determinants are positive and resolved.

Thus the finite-source data support two independent mean-response
directions inside the space of perturbations invisible to instantaneous
rank/Euler summaries. A single **fixed** response vector multiplied by
an arbitrary scalar control strength would have det J=0, unlike the
observed pattern. This is not a count of continuum fields and does not
exclude a general prefix-dependent scalar latent variable, whose averaging
can itself produce a full-rank response.

The distinction matters here: within cell00 alone the canonical A
determinants are `6.53e-10 +/- 5.55e-10` and
`1.52e-9 +/- 6.45e-10`, weaker than the complete-population result.
Combining01 and10 can yield rank two even though only one source is
active within either individual prefix. The demonstrated rank is that
of the ensemble mean-response map, not two simultaneously resolved
directions at every prefix.

This exploratory rank calculation was proposed after seeing the fixed
shared-label response matrix; it is not a preregistered decision test.
The four determinant readouts are correlated, not four independent
confirmations. No high-dimensional inverse, fitted response model, new
sampling or independent-source claim is involved.

## Source and reproducibility

### Full thermal response: cancellation hidden by the lifetime integral

The [complete common-label curves](https://github.com/LightChainr/Matching-One/blob/dcd63ace69eebbe80591e332a26830ca70560a85/notes/p334-paired-euler-thermal-response.md)
are now available from these same signed histograms. The minus -> D A peaks
are `+3.05999e-4 +/- 7.17445e-5` at p=.63313 and
`+3.14221e-4 +/- 4.54777e-5` at p=.62934. Plus -> D remains weak across
its principal thermal features; a physical root is not inferred for it.

N425 minus -> D E has an earlier negative feature,
`-6.31893e-5 +/- 2.09007e-5` at p=.60629, and a later positive feature,
`+4.76565e-5 +/- 2.02398e-5` at p=.66784. Both paired birth responses
are positive at these extrema, and which one dominates changes. These are
paired contrast signs, not acceleration claims for both individual clocks.
The corresponding main lobe areas `-3.08404e-6` and `+2.60091e-6` nearly
cancel; the complete E integral is only
`-5.25797e-7 +/- 2.21120e-6`. A small lifetime-integral response can
therefore conceal a finite thermal redistribution in this data block.

The numerical exchange at .63793 +/- .00871 is descriptive, with pointwise
peak uncertainties rather than selection-adjusted confidence statements.
N325's earlier negative feature is only1.34 pointwise SE, so the same
resolved two-lobe structure is not claimed there. The complete curves,
root descriptors and areas share the original20 batches.

The [final joint response analysis](https://github.com/LightChainr/Matching-One/blob/b582015e64e2d8a59e591c4822b14dedaea58b0f/notes/p334-common-label-tangent-joint.md)
retains the old trigger/contact factors, this common-policy response,
explicit cross-geometry channels and the supplied rank LOO coordinates
under one coordinator. The complete integer histograms separately preserve
the covariance of any additional full-p readout; numerical peaks are not
claimed to be extra independent evidence.

### Immutable inputs

Definition/code ffb70969; complete new response and signed histograms
4db356e1b026853468f94d59d938895a2367ceb7. Both reuse original e32a8593
suffixes and 959a7fa2 contact coordinates on the same40k paired prefixes.
`scripts/p334_common_label_response_rank.py` consumes only the saved20
batch vectors, and stores all linear Jacobians plus nonlinear determinant
and angle leave-one-batch values in
`results/p334-common-label-response-rank/score.json`. The designated
covariance coordinator can append these LOO columns to the existing
same-source factors without pretending they are independent observations.

The [exact common-policy theory](https://github.com/LightChainr/Matching-One/blob/cfaae36c8a73d85947547021f5c8c7a271c4367c/notes/p334-common-label-euler-tangent.md)
adds cross-channel support rules and a source-null theorem: a null direction
of the joint-class mark covariance is a class-constant source, so its
finite policy remains uniform and every future response vanishes. It also
distinguishes raw orientation parity from the normalized D coordinate,
whose denominator changes sign under exchange.

### Coordinated next analysis, not yet a result

The parallel research task `01a05696-fed1-7493-bd19-d581c88a3996` has
reported assigning HZsCM6 to a genuinely different question: use distinct
quartets of the same00 prefix to estimate local-Jacobian determinant and
determinant-square U-statistics. That distinguishes within-prefix response
dimension from mixture-generated ensemble rank. Its result is pending and
is not inferred from this note. This team will not duplicate that assigned
calculation. Finite-strength common policies and time-dependent response
kernels remain possible subsequent directions, not completed measurements.
