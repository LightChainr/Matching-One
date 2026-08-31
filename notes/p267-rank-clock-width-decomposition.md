# The broadening is in the rank-clock response, not binomial smoothing

**Removing the canonical binomial readout leaves 96.895% +/- 0.507% of the
observed N100-to-N400 increase in the centered thermal variance.** The
microcanonical rank-step variance increases by `0.662166 +/- 0.076535` in
`z=N^(3/8)(p-p_ref)`. This is a positive shape change in the underlying
finite rank response, not an artifact of a wider observation kernel.

No new data are generated. Sources are N100
[`7b30648`](https://github.com/LightChainr/Matching-One/commit/7b30648be558df0652a7ff22143cc87ed399d042)
and the independent N400
[`3e01b49`](https://github.com/LightChainr/Matching-One/commit/3e01b495b5b637b0070705e37b4137a9a0ef0d8b).
The observable is fixed to `D_A=P4[A_top](4i)-P4[A_top](2i)`.

## Exact finite-N relation

Let `c_k` be its signed threshold coefficients from K1/K2, with
`sum c_k=0`, and put `f_j=sum_(k<=j)c_k`. Then `f_N=0`. The two profiles are

\[
S(p)=\sum_k c_k H(p-k/N)
     =f_{\lfloor Np\rfloor},\qquad
C(p)=\sum_k c_k\Pr[\operatorname{Bin}(N,p)\ge k]
     =\sum_{j=0}^Nf_j B_{N,j}(p).
\]

S is the **rank-step profile**, not a guessed deconvolution. Its value is
already present in the archive at every integer occupation rank. C is its
ordinary Bernstein/binomial canonical readout.

Write `W=sum_j f_j`, normalized signed weights `w_j=f_j/W`, weighted mean
`mu_J=sum w_j j`, and centered second moment
`V_J=sum w_j(j-mu_J)^2`. Exact polynomial integration gives

\[
\int C=W/(N+1),\quad\int S=W/N,
\]

\[
\mu_C=\frac{\mu_J+1}{N+2},\qquad
\mu_S=\frac{\mu_J+1/2}{N},
\]

\[
\boxed{V_C=\frac{V_J}{(N+2)^2}
+\frac{\sum_j w_j(j+1)(N-j+1)}{(N+2)^2(N+3)}},
\qquad
\boxed{V_S=\frac{V_J+1/12}{N^2}}.
\tag{1}
\]

The `1/12` term is exactly the within-bin uniform second moment of the
piecewise constant step profile. The canonical formula follows from the
normalized Bernstein kernel `Beta(j+1,N-j+1)`: its conditional mean is
`(j+1)/(N+2)` and its variance is the second term in (1). This profile-kernel
Beta law is distinct from the threshold-CDF identity
`Pr[Bin(N,p)>=k]=CDF_Beta(k,N+1-k)(p)`.

If all w_j are nonnegative these are literal mixture-moment formulas. For
the signed geometry contrast they remain algebraic identities without
requiring nonnegative weights. They do not turn f_j into the probability
law of K1 or K2, and signed second moments need not generally be positive.

In particular, **canonical minus rank-step is not just extra binomial
variance**. Canonicalization both contracts the rank-grid coordinate by
`N/(N+2)` and contributes a conditional-Beta term. Subtracting a generic
`p(1-p)/N` noise variance would miss the normalization and mean change.
Multiplying (1) by `N^(3/4)` gives centered z variances; centering removes
the arbitrary choice of p_ref.

## What the existing data say

| Centered z variance/component | N100 | N400 |
|---|---:|---:|
| canonical | 1.4399945 +/- 0.0156460 | 2.1233811 +/- 0.0738206 |
| **rank-step, before canonicalization** | **1.4325846 +/- 0.0164411** | **2.0947509 +/- 0.0747484** |
| net canonical-minus-step correction | 0.0074099 +/- 0.0007955 | 0.0286302 +/- 0.0009353 |
| conditional-Beta term | 0.0632922 +/- 0.0001589 | 0.0494678 +/- 0.0002199 |
| contracted rank-grid contribution | 1.3767023 +/- 0.0158027 | 2.0739133 +/- 0.0740065 |

The rank-step increase is `0.6621662 +/- 0.0765352`, about 8.65 SE.
Of the canonical increase, the fraction retained in the rank-step response
is `0.9689482 +/- 0.0050742`; the net canonicalization correction contributes
only `0.0310518 +/- 0.0050742`. These fractions retain their strong paired
covariance, rather than dividing independent error bars.

The actual conditional-Beta term **decreases**, while the deterministic
contraction of the rank grid also weakens. Their net difference explains
why the small canonical-minus-step correction increases despite a narrower
conditional observation kernel. The primary broadening remains in f_j.

Each scale's 200/400 common-batch LOO calculations propagate the signed
area normalization, rank mean, and centered moment together. Cross-scale
covariance is block-diagonal because the counter domains are independent.
Every canonical/rank/component comparison within one scale is paired.

## A new working fingerprint, explicitly selected after these two scales

In the alternative coordinate `x=N^(1/4)(p-p_ref)`, the rank-step centered
variances are

```text
N100: 0.4530230 +/- 0.0051991
N400: 0.4684005 +/- 0.0167142.
```

Equivalently the two-area effective width is

\[
\gamma_{\rm eff}=\frac38-
\frac{\log(V_{z,400}/V_{z,100})}{2\log4}
=0.2379604\pm0.0135194.
\]

This is a **finite-regime, signed-profile width fingerprint**, not a new
critical exponent. The quarter-power coordinate was considered after both
scales were seen. Two centered moments do not show a full profile collapse;
means and other shape coordinates can still move.

It does give a concrete next-scale discriminator. Anchoring both models
to the same N400 rank-step variance, at N900:

| Conditional width hypothesis | Predicted N900 centered z variance |
|---|---:|
| p width proportional to N^(-1/4) | 2.5655354 +/- 0.0915477 |
| fixed profile width in z=N^(3/8)(p-p_ref) | 2.0947509 +/- 0.0747484 |

The quoted uncertainties are shared-anchor prediction uncertainties,
not target sampling errors, and are perfectly correlated between these
two fixed-ratio predictions. N900 has not been sampled in this analysis.
This is a conditional working target to freeze before future production,
not evidence already supporting either full profile mechanism.

## Artifacts and science card

Run `python3 scripts/p267_rank_clock_width.py`; output is
`results/p267-rank-clock-width/{score.json,REPORT.md}`. The JSON saves source
hashes, both full within-scale covariance matrices, all LOO vectors,
paired fractions, and the working N900 prediction covariance. A small
exact signed-weight test checks (1), including a negative weight.

- **Changed mechanism space:** canonical smoothing is quantitatively too
  small to explain the broadening; the finite rank-clock response carries
  the change. This is a physical decomposition, not another model p-value.
- **Not established:** no new critical exponent, no full profile collapse,
  no additional field count, and no assumption that the signed profile is
  a probability distribution.
- **Observer / geometry:** D_A ordinary P4 rank response, homothetic
  `2i,4i,1/2+i` families at N100/N400; no geometry re-selection.
- **Dependency:** the same original two independent archive blocks; this
  decomposition does not create independent evidence from reused moments.
- **Next discriminant:** independent-scale centered rank-step width, with
  the quarter-power versus fixed-critical-width predictions named above.
