# Occupancy mixing and rank selection in the global source response

## Technical summary

The original global source response now has a measured, additive decomposition.
For bulk `s=C_B+C_W`, **occupancy mixing contributes +.02432 to +.02478 to
the matching-root derivative; rank selection at fixed occupancy contributes
another +.00409 to +.00458** across the six finite norm-4 tori. Both components
are resolved. After following the moving root, both reduce rank-one population:
occupancy contributes −.06377 to −.06761, and rank selection −.02596 to −.03334.

The original H4 source derivative remains unresolved. At the two million-mark
endpoints its large central value and uncertainty sit mainly in the rank-selective
component, while the occupancy component is small and unresolved. Thus strong
conditional winding response, resolved rank-population response and unresolved
global H4 response are three distinct observations.

This calculation took **.859929 seconds**, with a150-dimensional full covariance,
three original source groups and100 aligned omissions per group. No configuration
replay, new counter, root finder, server action or research test suite was used.
[The JSON](latest.json) contains every component, chain contrast and covariance.

## Which source is visible to the global observer?

At each geometry g and occupancy K, let `r=q+1`, `E=q²`, and define under the
unperturbed microcanonical measure

```text
m_g(K) = E[s | K,g],
s_rank = E[s | K,r,g] − m_g(K),
s_spatial = s − E[s | K,r,g].
s = m_g(K) + s_rank + s_spatial.
```

For any function of `(K,q)`, conditional expectation makes the response to
`s_spatial` identically zero. This holds for all Bernoulli parameters p, hence
also for thermal derivatives, matching-root motion and the original global
`U=N^(13/8) P4[E_p]/(2D)`. Its source derivative therefore obeys

```text
v_total = v_occupancy + v_rank,     v_spatial = 0.
```

This extends the specific rank1-centered source result in
[the completed PR509 angular bridge](https://github.com/LightChainr/Matching-One/blob/fb01c44aa45e4f8d37d52144e2ad7c4adfe6ce40/experiments/p154-spatial-localization-20260831/REPORT.md).
The new information here is the actual nonzero allocation between the two
globally visible components. The residual zero is an identity, not another
experimental null result or a claim that spatial configurations are irrelevant.

## Root motion and topology redistribution are separately resolved

At N260, the root derivative decomposes as
`+.02473491±.00003104` (occupancy) plus
`+.00419573±.00006955` (rank), giving the original
`+.02893064±.00008134` per unit bulk source strength.
At N340 it is `+.02478400±.00003147` plus
`+.00409098±.00007119`, giving `+.02887498±.00008526`.
The rank contribution is resolved at least20.38 SE at every one of the six sizes.

For the root-comoving rank-one probability at N260, the corresponding derivatives
are `−.06532907±.00009603` and `−.03249425±.00062452`, giving
`−.09782332±.00061164`. At N340 they are
`−.06376945±.00009314` and `−.03334090±.00067363`, giving
`−.09711036±.00068302`. The rank-selective population term is resolved
at least21.49 SE across all six sizes.

These components share covariance; their errors must not be added independently.
The occupancy projection is the generally nonlinear geometry-specific function
`m_g(K)`, not a common affine-K thermal clock. A nonzero occupancy contribution
does not contradict the earlier rejection of that much narrower clock model.

## H4 uncertainty remains in a different readout

At N260 the original H4 source derivative is
`v=−5.97229±4.27362`. Its decomposition is
`v_occupancy=+.36379±1.30189` and `v_rank=−6.33608±4.34636`.
At N340 it is `v=+11.85937±9.19815`, from
`v_occupancy=+.91698±1.92619` and `v_rank=+10.94240±8.49261`.
Neither component supplies a new resolved endpoint H4 response. Do not turn
these weak totals into percentage attribution or sum the two components as
independent detections.

The frozen q2 `(1,−3,2)` and Jordan `(1,−2,1)` source-derivative contrasts were
also decomposed with full covariance along65→130→260 and85→170→340.
All component contrasts remain below one SE in absolute value. They do not
select a new source-rigid model, reverse the old high-statistics unperturbed
q2 rejection, or justify a fitted exponent.

## Exact calculation from existing moments

Write each fixed-K profile as
`x=<q>, e=<E>, m=<s>, tq=<qs>, tE=<Es>` and let `w_K=Bin(N,p)`.
The two visible q/source covariances are

```text
Jq_occupancy = Σ w_K x_K m_K − (Σ w_K x_K)(Σ w_K m_K),
Jq_rank      = Σ w_K (tq_K − x_K m_K).
```

Replacing x/tq by e/tE gives the two E/source covariances. These formulas avoid
division by rare rank probabilities. Thermal differentiation changes only w;
the microcanonical projection is fixed while differentiating p.

For any source component X, set `B=P4[E_p]`, `H=P4[E_pp]`,
`T=mean_g q_pp`, `D=mean_g q_p`, and `A=N^(13/8)/2`. Then

```text
p0dot_X = −mean_g Jq_X / D,
v_X = A [ (P4[(JE_X)_p] + p0dot_X H)/D
          − B (mean_g (Jq_X)_p + p0dot_X T)/D² ],
rank1dot_X = −mean_g JE_X − p0dot_X mean_g E_p.
```

All four v terms are saved separately. These are bulk derivatives, with no
additional factor N. The total and both projected sources use the same saved
pooled p0 in each central/omitted sample.

For interpretation only, a three-state representation
`E[s|K,r,g]=a(K,g)+b(K,g)q+c(K,g)E` exists on supported states. Its
coefficients can vary arbitrarily with K, geometry and size. This representation
does **not** identify a constant q/E field or a continuum operator. Empty
empirical states require no imputation in the formulas used here.

## Dependence, scope and reproducibility

The four cyclic sizes use100k old source-marked permutations each; N260/N340
use1M each. Every K row is a prefix of an existing permutation. The four cyclic
sizes share one random group; the two endpoints each have their own. Their
endpoint batches remain the original1000-counter block union its9000-counter
increment, with paired geometry deletion. No independent complement is used.
Each jackknife omission re-estimates `m_g(K)` and both projected responses, so
projection-estimation uncertainty is propagated rather than treating a fitted
projection as an independently fixed intervention.

This is a first-order, empirical source decomposition. Exponentiating a
conditional source mean generally does not reproduce the original finite-strength
measure. No causal birth-versus-duration allocation, asymptotic exponent,
cross-microscopic field identity or new independent sample is claimed.

Input/result identity is recorded in
[the contract](../../analysis/norm4_global_source_projection_contract.json).
Code was committed at `8b0aa3e4f507977ff6c3a98f0d1d03f0c66d65fc` before execution;
[the script](../../scripts/analyze_norm4_global_source_projection.py) records
all input and code SHA256 values and the actual Python environment. Across the
central estimate and all aligned omissions, component addback error is at most
`5.21e−11`; the old global source derivative is reconstructed within `6.59e−11`.
These are calculation identities, not extra scientific evidence.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
python scripts/analyze_norm4_global_source_projection.py
```

The requested primary artifact is this repository Markdown/JSON pair. No extra
HTML dashboard or chart is needed: the exact projection equations and two
endpoint decompositions carry the distinction, with the complete six-size
numerical record retained in JSON.

## Next mechanism-changing move

Prioritize a **specified source acting before a later rank observation**, and
predict its subsequent rank-population and original-U response together. The
completed P334 finite common-policy experiment already demonstrates the logic:
instantaneous rank/Euler distribution can be preserved while future birth
responds. That dynamic statement and the static conditional-projection zero
are compatible, because they concern different observation times.

For norm4, the missing bridge is a paired, lag-resolved source-to-birth/exit
kernel with an explicit source time and later readout time, or a concrete
microcanonical mechanism predicting the measured fixed-K rank covariance.
Same-time marginal source/line profiles do not contain that joint kernel.
Reuse an existing two-time archive if it has the required marks; do not relabel
the present static projection as that experiment. Stronger centered line
association alone cannot identify the original global H4 source.
