# P334: exact birth-age memory and a separate direct-birth channel

Date: 2026-08-30. Status: **exact finite-volume results; conditional scaling conjecture kept separate**.

This review reads main at `4d70c1787ff97dbb98cb5e96022f947bb8fad97e` and follows #334 through `2b5d706d2e91fd0d4119416c2f9ff82708ad2297`. Projective birth marks, source/sink fluxes, constant integral saturation, the BA/TM hazard program, and corrected translation/relative-phase Hall compression already exist. They are not proposed here as new. This note does not solve the remaining arbitrary-HNF Hall inequality.

## 1. The two-time rank source determines the unmarked birth-pair law

Under common independent uniform site priorities, let `tau1<=tau2` be the first and second essential rank births, allowing equality, and `Fi(p)=P(taui<=p)`. Pathwise,

```text
r(p)=1[tau1<=p]+1[tau2<=p].
```

For `p<=q`, multiplication gives the exact identity

```text
H(p,q):=P(tau1<=p,tau2<=q)
       =E[r(p)r(q)]-F1(p)-2F2(p).
```

For `p>q`, ordering gives `H(p,q)=F2(q)`. Thus the full two-time rank correlation, together with marginal curves, determines the entire unmarked birth-pair distribution, including its diagonal mass.

This gives an intrinsic extension of #337:

```text
Z(p,q;s,t)=E[exp(s(r(p)-1)+t(r(q)-1))].
```

The mixed source derivative recovers the two-time moment after adding known marginal terms. By contrast, the one-time source algebra and arbitrary p derivatives contain only one-time information. The current sparse `(K1,ell,K2,DIRECT_RANK2)` archive already supports the relevant joint analysis. These are correlated reanalyses, not independent new samples.

## 2. Actual history dependence: an exact 1/57 witness

Take nearest-neighbor square-site percolation on the Gaussian quotient `3+i`, with period columns `(3,1),(-1,3)` and `N=10`. Grow a uniformly random site permutation. At `k=5`, condition on ambient rank one and primitive projective line `ell=(1,0)`.

| K1 | ordered prefixes | exit-weighted prefixes | P(K2=6 given k=5,ell,K1,K2>5) |
|---|---:|---:|---:|
| 4 | 1440 | 1920 | 1920/(5*1440)=4/15 |
| 5 | 4560 | 6480 | 6480/(5*4560)=27/95 |

The difference is exactly `27/95-4/15=1/57`.

Therefore `(current occupied count, rank, projective line)` is not a Markov state for the observed process under the actual uniform-permutation law. This is not merely a failure of strong lumpability between two microscopic states.

### Counting certificate

For rank-one `S`, `|S|=k`, let `a_j(S)` count ordered prefixes ending at S whose first birth is j. Delete the last vertex:

```text
a_j(S)=sum over v in S of:
  (k-1)! 1[j=k],   if r(S-v)=0;
  a_j(S-v),       if r(S-v)=1.
```

A rank-two predecessor is impossible. Every rank-one predecessor has the same line. The certificate checks `sum_j a_j(S)=k!` for every rank-one subset.

With `x(S)` the number of vacant vertices whose insertion produces rank two,

```text
h(k,j,ell)=sum_S a_j(S)x(S)/((N-k)sum_S a_j(S)).
```

All **30,240 ordered length-five prefixes** were also enumerated independently of the recurrence, reproducing all four integers in the table.

### This does not contradict the existing uniform-layer hazard results

There are 50 rank-one states on this line at k=5. Their uniform exit hazard is `xi_5=7/25`. The two history classes have weights `6/25` and `19/25`; their conditional hazards average back to `7/25` exactly.

More generally, under uniform fixed-line layer measure mu, set `w_j(S)=a_j(S)/k!` and `h_x(S)=x(S)/(N-k)`. Then

```text
h(k,j,ell)=xi_k+Cov_mu(w_j,h_x)/E_mu[w_j].
```

This history-conditioned Radon--Nikodym identity is compatible with the repository's uniform-layer monotonicity, BA/TM conditions and ULC conjecture. It is not a counterexample to any of them.

### Exact age augmentation, not an invented latent mode

While rank is one, its observed rank-line history is completely specified by entry time K1 and the unchanged line ell. Hence `(k,r,ell,K1)` is an exact time-inhomogeneous Markov representation of this **observed** history. Suppressing K1 gives an age-dependent survival process. These coordinates do not determine the microscopic configuration.

In continuous priority coordinates, define

```text
rho_ell(t,p)dt=P(tau1 in dt,tau2>p,plateau line=ell), t<p;
d_p rho_ell(t,p)=-h_ell(p|t)rho_ell(t,p).
```

The direct `0->2` flux is a separate line-free channel. Whether birth-age dependence survives near-critical scaling remains unproved here. A deterministic increasing intrinsic clock such as `u(p)=(F1(p)+F2(p))/2` can remove clock conventions, but cannot remove genuine history dependence. A data-estimated clock must retain joint uncertainty.

## 3. Direct births measure pathwise quadratic variation

Define the unnormalised directed subset-edge count

```text
E02[k]=#{(S,v): |S|=k, v not in S, r(S)=0, r(S+v)=2}.
```

Then the exact Bernoulli-priority flux and mass are

```text
j02(p)=sum_k E02[k]p^k(1-p)^(N-1-k),
D_N=P(K1=K2)=integral_0^1 j02(p)dp
   =sum_k E02[k]/(N*binom(N-1,k)).
```

Dividing the raw total by the number of directed subset edges would use the wrong path measure.

| Gaussian generator | N | direct subset edges | D_N |
|---|---:|---:|---:|
| 2+i | 5 | 0 | 0 |
| 3 | 9 | 45 | 3/35 |
| 3+i | 10 | 80 | 5/63 |
| 3+2i | 13 | 793 | 304/3465 |
| 4 | 16 | 4624 | 2809/45045 |
| 4+i | 17 | 8823 | 3511/60060 |

The N17 count 8823 and rank-one axis/diagonal counts 36516/2380 match the existing `e34140d` census. This is an explicit cross-check, not a revalidation of the entire repository. No asymptotic fit is made from these tiny geometries.

Every path has either two unit rank increments or one double increment. Consequently,

```text
sum_k(r_(k+1)-r_k)^2=2+2*1[K1=K2].
```

For fixed finite N and interior p,

```text
lim_(h down to 0) E[(r(p+h)-r(p))^2]/h=M'(p)+2*j02(p).
```

Two distinct site priorities fall in a length-h interval with probability O(h^2), proving the limit. Pathwise quadratic variation is not the same observable as the one-time source moment `E[(r(p)-1)^2]`.

All six censuses satisfy `int j01=int j12=1-D_N` and `int M'=2` exactly.

## 4. Conditional six-arm prediction: NOT a proved result

A useful geometric lemma to attempt is that a direct rank-zero to rank-two birth on a nondegenerating honest torus requires three occupied arms and three separating dual vacant arms to a fixed fraction of the injectivity radius. A matching lower-bound gluing construction is also needed. Neither is proved here.

Assume additionally per-site comparability to the six-arm probability within the near-critical window, a nonzero scaling amplitude, and adequate control of tails outside the window. Using the triangular-site exponent input gives

```text
alpha6(L)=L^(-35/12+o(1)),
window width=L^(-3/4+o(1)),
j02(pc)=L^2 alpha6(L)=L^(-11/12+o(1)),
D_N=L^(-5/3+o(1))=N^(-5/6+o(1)).
```

A fixed-ratio limit `D_(4N)/D_N -> 2^(-5/3)` requires the stronger regular-scaling hypothesis, not only an exponent with o(1). Importing the triangular exponent to square-site percolation is another universality hypothesis. The direct-birth/arm identification itself remains a separate gap.

Even a suitable upper bound could exclude simultaneous essential births in the scaling limit. None of this identifies the matching H4 residual with the collision channel or supplants the Q4/Jordan candidate.

## 5. Literature and next scientific discriminator

Camia--Feng, arXiv:2508.16047v2 (revised 2026-06-01), provides triangular-site lattice logarithmic-pair two- and three-point limits. Roux--Ribault--Jacobsen, arXiv:2604.24491 (2026-04-27), relates torus one-point functions to sphere four-point functions at another central charge. Ang et al., arXiv:2604.05503 (2026-04-07), proposes exact loop-model three-point formulas with multiple checks. Camia--Foit--Nivesvivat, arXiv:2605.04395 (2026-05-06), supplies new FK pivotal spatial-density formulas in a bulk-boundary CFT/SLE setting.

These results motivate jointly comparing one-point, two-time and explicitly identified local/defect insertions. None supplies the missing lattice-observable dictionary merely from a fitted logarithm or modular covariance.

The proposed decision test is to separate **one-time topology probabilities, birth-age dependence, and direct-birth diagonal mass**. If the last two vanish after scaling, a memoryless projective limit becomes plausible. If collisions vanish but age dependence persists, a simple but memoryful birth process is the natural target. Finite-volume path non-Markovianity does not refute a finite-dimensional spatial RG or Jordan representation.

## Reproduction and scope

Run:

```sh
python scripts/p334_birth_age_collision_review_20260830.py --output /tmp/p334-birth-age-collision.json
```

The standard-library script enumerates all subsets of the six specified Gaussian quotients, checks the exact masses and N17 reference counts, and independently enumerates 10P5 prefixes to verify the memory witness. It was executed locally. Separately implemented Python/C++ checks and six local regression tests also passed in the accompanying analysis bundle. Repository-wide CI and production simulations were not run. The script is a bounded verification oracle, not a new production stack.

References:
- https://github.com/LightChainr/Matching-One/issues/334
- https://github.com/LightChainr/Matching-One/commit/e34140d43b5999bb805cbb1735f3ccc2b2ed5c4c
- https://github.com/LightChainr/Matching-One/commit/874142c649afc686d1b5839a18bb38f5bb0486b5
- https://github.com/LightChainr/Matching-One/commit/2b5d706d2e91fd0d4119416c2f9ff82708ad2297
- https://arxiv.org/abs/2508.16047
- https://arxiv.org/abs/2604.24491
- https://arxiv.org/abs/2604.05503
- https://arxiv.org/abs/2605.04395
- https://arxiv.org/html/1305.5526 (triangular six-arm estimate (6.10), four-arm exponent, near-critical length scaling)
