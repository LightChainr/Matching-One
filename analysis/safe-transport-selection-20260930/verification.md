# Narrow verification: safe transport and survival selection

2026-09-30. **No actionable correctness error found in the reviewed source.**
The random-order rejection/first-proposal coupling, endpoint eligibility,
common-cell decomposition and aligned batch covariance implement the stated contract.
This review supplies no new precision estimate and does not certify a cloud run.

Reviewed at source commit `5054c02383242473834045fcbddc7379d83acfb9`;
the six reviewed files were unchanged through observed HEAD
`a58643d1f787cd73e3f12610d517b2fc5d0e3b91`.
Only this verification file was written; no commits or pushes.

## Algebra and source review

| Item | Finding |
|---|---|
| First-birth source | engine.cpp:18–27 generates the original complete Fisher–Yates permutation from its archived seed, then stops at the first positive rank or b. It includes every rank-one entrance J1≤b, including prefixes whose **old** endpoint later reached rank two. Direct 0→2 entrances and no entrance by b are separately recorded. Restricting to old rank_b=1 would have selected the source twice; the code does not do that. |
| Separate continuation stream | The old permutation RNG and new `future` RNG are separate (engine.cpp:14). run.py derives one continuation seed per original batch from the exact contracted SHA256 namespace, first eight bytes in big-endian order. Conditional branches share their original batch dependency; they are not new independent preparations. Birth-only alignment to the old CSV is appropriate: the new endpoints are deliberately different realizations. |
| Random-order rejection | At descending position r, swapping with a uniform index in [0,r] draws uniformly from the as-yet unexamined vacancies. Previously rejected sites remain in the frozen suffix for this step. On acceptance, swapping that site with `back()` and popping removes exactly the accepted site, preserving every rejected vacancy. Starting the next partial shuffle from any resulting array again gives a uniform random order; its arbitrary ordering cannot bias the next accepted site. |
| Local completion query | Geometry.completes_if_inserted:94–107 uses the same lifted-neighbour potential conflict as completion_count:70–90. It checks neighbours in a common old component, detects a winding transverse to the current rank-one direction, and never inserts a site. `find` may compress union-find paths; this is a representation change, not a change in occupation or homology. The query is called only on current vacancies of a rank-one geometry. |
| First-proposal flag | For each safe v, P(accepted v)=1/s and P(first proposal safe, accepted v)=1/m. Consequently the first-safe indicator has probability s/m independently of accepted v, conditional on the geometry. engine.cpp:36 updates the permanent flag only on the first proposal; later rejected proposals are not additional natural insertion times. If the flag was already zero it remains zero. `first_rejection_count` records the first killing **insertion count** k+1, rather than a tally of all rejections. |
| Whole-path coupling | The above conditional probabilities hold for any vacancy-array ordering and past flag history. Thus, conditional on the accepted safe path, E[I_end]=∏s_i/m_i=W. Alive accepted transitions have probability 1/m, exactly the killed original uniform kernel; dead flags do not stop the reference chain. No completion-census estimator or product of estimated probabilities enters this construction. |
| Cemetery | If s=0, all m vacancies are examined once, no site is inserted, the code sets rank=-1 and alive=0, and exits the count loop. There is no retry cap or infinite rejection loop. This branch was source-reviewed, not separately simulated. |
| Endpoint and direction | Successful paths contain exactly b distinct occupied sites. Safe insertions preserve the primitive direction; the explicit endpoint direction assertion justifies keying endpoint cells by `dx_entry,dy_entry`. The pair counter is invoked once at b and supplies c and e; Y=2e/(N−b−c). No entrance, direct rank-two entrance, reference cemetery and c=N−b are counted separately, with c=N−b excluded because Y is undefined there. |
| Fixed parameters | run.py passes L=512 and b=155385; score.py uses a=154646, b=155385 and m=106759. These match contract.json. Source replay, cohort definitions and endpoint eligibility refer to the same fixed count clock. |

For each eligible cell C=(D,c), write nS, sumY, nN=sum I and sumIY
for the reference count, reference outcome sum, natural-alive count and
natural-alive outcome sum. score.py:95–101 accumulates precisely these
four quantities separately for each cohort. In that cell,

```text
mu_safe = sumY/nS
mu_natural = sumIY/nN
Cov_emp(Y,I)/mean_emp(I)
  = [sumIY/nS - (sumY/nS)(nN/nS)] / (nN/nS)
  = mu_natural - mu_safe.
```

The covariance here uses the empirical probability convention (denominator nS),
not an nS−1 unbiased sample-covariance convention. The scorer implements the
difference directly, so no erroneous finite-sample covariance factor is introduced.

score.py:31–43 retains only cells with nN_early>0 and nN_late>0;
these automatically also have positive reference counts. It uses the **same**
weight w=nN_early*nN_late/(nN_early+nN_late), support and normalization
for natural and reference contrasts, then sets selection=natural−reference.
Thus the decomposition is of a natural-overlap-weighted endpoint-cell contrast,
not an unconditional population contrast or a unique causal attribution.
Reference observations with dead flags remain in their reference cell means.
Natural conditional means have the original survivor law in distribution;
they need not reproduce the previously observed archive endpoints row by row.

On each deletion, evaluate() reconstructs cell aggregates from the remaining
eligible batch rows, redetermines natural-cohort overlap, and recomputes all
means and weights. Precomputed per-row eligibility is fixed by the endpoint
contract, so excluding the deleted batch's eligible rows is the same operation
as reevaluating that fixed eligibility rule. If no natural overlap remains,
the code raises `not_scoreable` rather than inventing a zero contrast.

score.py:107–110 uses the 14 aligned three-vectors and
`(13/14)*sum((theta_minus_batch - deletion_mean) outer itself)`.
The full 3×3 covariance, including cross terms, is correct; its singularity
is expected because the third coordinate is the first minus the second on
every deletion. The scorer saves all 14 deletion vectors. It does not add
component variances as if the natural/reference estimators were independent.

## Actual check 1: one deterministic L4 vacancy-query control

Used one fixed square-L4 geometry, v=x+4y:

```text
A={0,1,2,4,5,8,12}; rank=1; D=(0,1); vacancies=9; c=1.
```

For each of its nine vacancies, the **actual current C++** Geometry was copied
and `copy.insert(v)==2` compared with `g.completes_if_inserted(v)` on the
original geometry. All **9/9 agreed**; only site 3 completes. The original
direction and completion_count remained unchanged after all queries.
This single geometry therefore exercises both safe and completing answers.

To keep the sole-file write boundary, the Geometry source was read into memory,
followed by a tiny C++ harness; `/usr/bin/clang++ -std=c++17 -O1 -S -emit-llvm
-x c++ -o - -` took stdin and returned LLVM IR on stdout. Research Python
3.11.15 and llvmlite/LLVM 20.1.8 executed it in memory; no test source or
build artifact was saved. The first JIT setup aborted before executing the
control because Apple Clang 21 emitted the unsupported backend attribute
`"probe-stack"="__chkstk_darwin"`. Removing only that attribute from the
in-memory IR allowed the single fixed control to execute successfully.
No Geometry logic, repository source or production binary was modified.
This is a local semantic control, not a verification of the cloud compiler.

## Actual check 2: one exact analytic coupling/decomposition control

Used a minimal two-step monotone abstract kernel with vacancies
`{a,b,c,d,x}`: x always completes and {a,b} is a synergy pair. Thus initially
m=5,s=4; after a or b there are two safe sites among four vacancies, and
after c or d there are three. This is an analytic kernel control, not an
extra lattice preparation or Monte Carlo replica.

A Python standard-library Fraction recursion enumerated the descending
partial Fisher–Yates choices, first-proposal flag update and swap/pop removal
literally as in the engine; no RNG was used. Every branch deleted only its
accepted site and retained all rejected vacancies. Exact results:

- Each initial safe site has P(accepted v)=1/4 and
  P(first-safe, accepted v)=1/5.
- For every complete accepted path, E[I_end | path] equals W:
  2/5 for paths starting a or b; 3/5 for paths starting c or d.
- With f=1{endpoint={c,d}}, E_safe[f]=1/6, E[I]=1/2,
  E[If]=1/10 and E[f | I=1]=1/5.
- `natural−safe = 1/5−1/6 = 1/30 = Cov(f,I)/E[I]` exactly.

No other physical configuration, lag grid, census, archive score or production
was run. Support and batch covariance were reviewed algebraically from the
current scorer; no synthetic batch-replica suite was added.

## Scope and handoff

Read the requested note, engine, contract, runner, scorer and added Geometry
query, with their directly relevant reuse interfaces. No SSH, cloud operation,
production, archive replay, extra precision calculation, commit or push.
No substantive issue requires delaying the main investigator's bounded job.
Cloud execution, runtime row correspondence, empirical support and numerical
results remain the main thread's responsibility; they were not claimed checked.

Reviewed SHA256 values:

```text
notes/safe-path-selection-coupling-20260930.md
4001b5127d4ceb7b362d6f772a000e9c8a4cccbab102004e0af133f82c42fc96
analysis/safe-transport-selection-20260930/engine.cpp
4d5d1ed48beffc0e13f8db12da67024118a4b1ea75fa7f8e68b3c6153c072c0b
analysis/safe-transport-selection-20260930/contract.json
22ca3b823989df75eeb1fe6091bdeee38eaf6229950c4cc441421fcbd3850bc2
analysis/safe-transport-selection-20260930/run.py
a5bf404c70bf65cd28dc18a4aefa2b59729f65c64b74c11b0f0f74a548e9f31c
analysis/safe-transport-selection-20260930/score.py
50e681e2d1422e0d2b1c4aa089bdc462d1a5b7767de2f585821c50a1db2d7b74
analysis/completion-hazard-production-20260929/engine.cpp
da45ae8a2682f9edb24eaf87a6733c55013f1fb95b3622bf2faf97208c405c0b
```
