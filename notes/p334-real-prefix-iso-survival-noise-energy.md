# P334: a flatter real-prefix clock can have stronger spatial noise energy

At the same survival level1/2, the original real N425 witness B has a
**smaller mean-response slope but30.6% larger first-degree noise energy** than
witness A. Its clock-only slope prefactor is27.1% lower; its pivotal
concentration is79.3% higher. The spatial concentration more than compensates
for the flatter clock.

This uses only the already-computed full clocks6358ba49 and all-site pivotal
polynomials1c06230b. No network DP, new suffix, new MonteCarlo sample, or
cohort expansion was performed. These are the original counters43042514269
and43042505280, not the later15-direct-gate pair.

## Compare equal response levels, keeping the occupation parameter typed

For each fixed rank-one prefix let d=173 be the number of remaining sites and
f_k the number of safe k-subsets. Under independent occupation u of these
remaining sites,

\[
S(u)=\sum_{k=0}^{d}f_k u^k(1-u)^{d-k}.
\]

Choose its unique median u_* with S(u_*)=1/2. This u is **remaining-site
occupation**, not the original full-N canonical p and not k0/N. Each fixed
prefix has its own u_*, so both indicators have variance1/4 at evaluation.

The saved D_v polynomial counts other-site subsets on which adding v causes
birth. Thus its fixed-u pivotal probability is

\[
I_v(u)=(1-u)^{d-1}D_v\!\left(\frac{u}{1-u}\right),\qquad
\chi(u)=-S'(u)=\sum_v I_v(u).
\]

In the orthonormal u-biased product basis, the singleton safe-indicator
coefficient is `-sqrt(u(1-u)) I_v(u)`, so

\[
E_1(u)=u(1-u)\sum_v I_v(u)^2
      =\underbrace{u(1-u)\chi(u)^2}_{\text{clock-only prefactor}}
       \underbrace{\sum_v\left(\frac{I_v(u)}{\chi(u)}\right)^2}
                   _{C_{\rm pivotal}(u)}.
\]

The whole first-order energy is therefore not just the squared clock slope.
The second factor records how that slope is distributed across actual sites.

## Two real-prefix readouts

Both source configurations are N425 second HNF[[425,268],[0,1]], k0=252,
age10, ell=(12,-19), with no original single-site direct gate. All values
below are deterministic high-precision evaluations of exact stored integer
polynomials; they are not MonteCarlo estimates with sampling error bars.

| Readout at S=1/2 | A:43042514269 | B:43042505280 |
| --- | ---: | ---: |
| Remaining-site median u_* | .09136112631 | .10638256970 |
| Clock slope χ | 6.70764011 | 5.35008304 |
| Clock-only prefactor u(1−u)χ² | 3.73501426 | 2.72109155 |
| Pivotal concentration C | .02671664 | .04789149 |
| Effective pivotal sites1/C | 37.42986 | 20.88054 |
| Positive pivotal support | 83 | 127 |
| **First-degree energy E1** | **.09978703** | **.13031713** |
| Total indicator variance | .25 | .25 |
| **Sum of degree≥2 energies** | **.15021297** | **.11968287** |
| E1 fraction of total variance | 39.9148% | 52.1269% |

B has more sites with nonzero influence, yet a smaller effective support:
its response is concentrated much more unevenly. The exact multiplicative
comparison, rounded for display, is

\[
\frac{E_{1,B}}{E_{1,A}}
=\frac{[u(1-u)\chi^2]_B}{[u(1-u)\chi^2]_A}
 \frac{C_B}{C_A}
=0.7285358\times1.7925716=1.3059526.
\]

Thus the change of energy ordering is specifically explained by the spatial
distribution of pivotality, not by a steeper unmarked response curve.

## Separate the uniform tangent from spatially nonuniform first-order energy

The normalized uniform singleton direction is
`d^(-1/2) Σ_v (X_v-u)/sqrt(u(1-u))`. Its projection energy is completely
determined by the clock:

\[
E_{1,\rm uniform}=\frac{u(1-u)\chi^2}{d},\qquad
E_{1,\rm spatial}=u(1-u)\sum_v\left(I_v-\frac{\chi}{d}\right)^2.
\]

These are orthogonal components, with E1=E1,uniform+E1,spatial. They are not
an empirical null model or a fitted regression decomposition.

| First-order component | A | B |
| --- | ---: | ---: |
| Clock-determined uniform projection | .02158968 | .01572885 |
| Spatially nonuniform excess | .07819735 | .11458828 |
| Spatial excess fraction of E1 | 78.3642% | 87.9303% |

The smaller clock projection in B would suggest weaker response if spatial
information were omitted; its larger nonuniform component reverses that
ordering. This is the missing information entering the first positive
noise-degree coefficient. For a fixed marginal u and replica correlation
rho, `Cov(f(X),f(Y))=rho E1+O(rho²)`: B has the larger initial noise-covariance
slope despite its flatter occupation-response slope.

Only E1 and the aggregate degree≥2 remainder are evaluated. No individual
E2,E3,... or P437 fifth-degree production amplitude is inferred. The separate
five-site isoclock counterexample established the general nonidentifiability
boundary; the present two real prefixes have different full clocks and serve
as an actual-source measurement of the slope/concentration split, not another
same-clock counterexample.

## Computation, dependency and handoff

`scripts/p334_real_prefix_median_noise_energy.py` evaluates F and D_v using
the existing research environment, mpmath80-digit arithmetic and180 bisection
steps. Root brackets and50-digit decimal outputs are saved in
`results/p334-real-prefix-median-noise-energy/iso_survival_median_energy.json`,
including each site's I_v and singleton energy. There are zero new network
solves. Numerical display precision is not population precision.

These two posthoc configurations share the original cooperative-closure
production dependency block. They are not a paired/stratum-weighted estimate
over the20k counter population, nor a replacement for that loading. The main
agent is separately performing the declared paired rank-one-stratum loading;
this task supplied source-contract/replay pointers and did not duplicate or
change that computation. No server was started or contacted.
