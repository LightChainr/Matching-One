## #437: same-energy fifth-difference representation, and a precise sparse-support obstacle

Follow-up to the 20k cost result `1751e1c`; **no new Monte Carlo**.
Branch `theory/p437-positive-difference-bridge-20260831` contains the compact
derivation, executable 32-vertex mixed-derivative kernel, and exact certificate.

1. With `a=9765/32768` and six dyadic nodes, the original six-level energy is
   exactly `a [t0,...,t5]K`. Hermite--Genocchi gives
   `A=a/5! E_Dirichlet K^(5)(T)`; hence
   `A=a E_T sum_|S|=5 <D_SF,T_T D_SF>`.
   All six coefficients and spectral multipliers through degree20 match exactly.
2. Low-degree removal is now **pointwise** in each fifth mixed derivative,
   eliminating the original million-fold low-order cancellation noise. But a
   single noisy pair product is still signed; only the conditional-mean
   square is pointwise positive. A degree6 exact counterexample returns -1.
3. **Do not launch uniform five-bond sampling.** There are
   `C(224,5)=4,493,032,544` candidates. On a single degree5 Fourier control,
   the uniform estimator has exact variance **399,010,376.139**, worse than
   the preceding six-level control. The known-support proposal has zero
   variance on the same control: this is a proposal problem, not a no-go for
   mixed differences.

**Next useful gate:** topology-aware candidate sets with exact inclusion
probabilities (full support or a proof of zero omitted derivatives), then a
small same-estimand variance measurement. Merely selecting currently pivotal
bonds is not sufficient: collective fifth differences need not contain a
single-bond pivotal at the starting configuration.

Observer/source/geometry unchanged: square-bond p=.5, N112 rho C3 Etop,
self-source spectral energy. This derivation is deterministic and does not add
an independent production evidence block or identify a physical field.
