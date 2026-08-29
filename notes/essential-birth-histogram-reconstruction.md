# Essential-birth views already present in the threshold-rank archive

Status: deterministic reinterpretation for Issue 269. No new samples are generated.

The committed axis `L=8` pilot contains 100,000 paired `(K_minus,K_plus)` observations, their two
marginals, integer moments, metadata, and checksums. The essential-birth filtration oracle gives these
historical endpoints a direct pathwise interpretation:

```text
K_minus = first ambient-H1 birth,
K_plus  = second ambient-H1 birth.
```

## Archive integrity

Before deriving anything, the analysis verifies every declared SHA-256 checksum, both marginal totals,
the joint total, joint-to-marginal reconstruction, every stored joint integer moment, and the support
condition `K_minus<=K_plus`. It also reproduces the archived `M` and `M_prime` evaluation independently.

## CDF and density reconstruction

For a random site permutation with continuous uniform priorities, let `tau_1=U_(K_minus)` and
`tau_2=U_(K_plus)`. Their empirical CDFs are obtained exactly from the threshold histograms by mixing
binomial upper tails. Write them as `F_1,F_2`, with densities `f_1,f_2`. Then

```text
P(R=0) = 1-F_1,
P(R=1) = F_1-F_2,
P(R=2) = F_2,
M       = F_1+F_2-1,
rho     = M'/2 = (f_1+f_2)/2.
```

Consequently the repository's canonical CDF `(1+M)/2` is the equal-weight mixture of the two essential
birth times. The finite matching root is the median of that mixture and the balance point `P(R=0)=P(R=2)`.

## Center and lifetime

Set `C=(tau_1+tau_2)/2` and `W=tau_2-tau_1`. Uniform order-statistic formulas convert the committed
joint integer moments into exact rational means and variances for `C` and `W`. In particular,

```text
integral_0^1 P(R=1) dp
  = E[W]
  = E[(K_plus-K_minus)/(N+1)].
```

Thus the historical neutral-area statistic is exactly the expected lifetime of the unique ambient-rank-one
phase, not merely an abstract gap covector.

## Irrecoverable marks

The joint endpoints are sufficient for the CDFs, densities, balance ratio, and center/lifetime moments.
They do not retain:

- the projective winding line `ell` during the rank-one plateau;
- the integral saturation index of the ambient subgroup;
- local marks at the first and second birth sites.

Those fields require a future production-compatible stream and cannot be reconstructed from the committed
endpoint histograms. This analysis reuses the same pilot data and is not independent evidence.
