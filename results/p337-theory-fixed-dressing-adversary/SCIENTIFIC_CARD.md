# Scientific card: fixed identity dressing is the strongest low-parameter adversary

- Mechanism changed: fixing both transfers to the theory values
  `lambda0=2^-13/8` and `lambda_id4=2^-21/8` and fitting only two amplitudes
  gives `q=0.803/2`, `p=0.669`.  It is not merely viable: its descriptive AIC
  `4.803` is lower than free-single `5.979`, free-lambda recurrence `6.077`, and
  rank-3 Jordan `6.084`.
- Exact dressing fingerprint: fitted amplitudes are `c0=-0.05528` and
  `c_id4=+0.04024`.  Because `lambda_id4/lambda0=1/2`, dressing/leading
  magnitude halves exactly each generation:
  `0.728 -> 0.364 -> 0.182 -> 0.091 -> 0.0455`.
- Rank-3 adversary: the minimal same-base form
  `lambda0^n(c0+c1 n+c2 n(n-1)/2)` also passes (`0.0839/1`, `p=0.772`).  It
  predicts N1360 `-0.0007896`, essentially identical to the free-lambda
  recurrence `-0.0007956`; those two mechanisms are forecast-degenerate.
- Other comparators: free-single remains acceptable (`1.979/2`, `p=0.372`),
  while fixed single H4 (`15.843/3`) and scale-neutral (`68.940/3`) fail.
- N1360 identifiability: full cross-model source covariance limits identity
  dressing versus rank-3, free-recurrence and free-single to maximum
  `0.594`, `0.695`, and `2.308 sigma`, respectively, even before adding new
  measurement noise.  A lone N1360 run therefore cannot resolve those classes.
- What N1360 can do: identity versus fixed-single has a source ceiling of
  `3.878 sigma`, but needs future measurement SE below `5.10e-5`; neutral is
  easy to reject.  This is not currently an efficient universal production.
- Observer/sector/source: same-lineage sign-aligned `K_A=d_eta log W_A` H4
  amplitudes, with full paired covariance inside each generation and independent
  covariance blocks across N85/N170/N340/N680.
- Does not prove: that the second contribution is literally the identity family
  rather than a Jordan or other correction state.  It establishes a sharper,
  parameter-economical theory adversary and the measurements needed to break
  the degeneracy.
- Next discriminator: add a geometry/modulus covector whose identity-family
  dressing and same-base Jordan response differ; simply increasing lineage N is
  source-covariance limited.
