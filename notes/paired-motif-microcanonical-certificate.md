# Paired motif microcanonical certificate

Status: exact-validation slice of Issue 40.

For any collection of `c_m` embeddings of an `r_m`-site motif on `N` labelled sites, a uniformly chosen
fixed-rank configuration with `K=k` occupied sites satisfies

```text
E[T_m | K=k] = c_m * binom(k,r_m) / binom(N,r_m).
```

Consequently, two same-N orientations with equal embedding multiplicity have an exactly zero-mean
difference at every rank. This does not depend on empirical recentering or a fitted coefficient.

The executable certificate freezes four families: NN edges, diagonal pairs, elementary faces, and
right-angle triples. It exhausts all configurations for the primitive Gaussian pairs `(2,1)/(1,2)` at
`N=5` and `(3,2)/(2,3)` at `N=13`. It also checks the expectation formula independently on axis `L=3`
and diamond `L=2`, compares direct and incremental last-activated-vertex counters, and verifies the full
joint motif-count histogram under a determinant-one period-basis change.

The two tiny conjugate/swapped pairs are configurationwise degenerate under the deterministic quotient
labelling, so they validate the counters and fixed-rank algebra but cannot demonstrate useful covariance.
The certificate therefore also checks equal multiplicities for all five declared production pairs
`N=65,85,130,145,170` and records a deterministic four-site witness on each pair where every frozen motif
difference is nonzero. Their fixed-K means are still exactly zero by the proved multiplicity formula;
only their covariance with the target remains empirical.

This result establishes admissible exact zero-mean controls only. It does not show that the controls are
correlated with the orientation target, reduce variance, improve variance per wall time, or pass the
Issue 40 production gate. Those questions require pilot-frozen or cross-fitted coefficients on fresh
paired samples.
