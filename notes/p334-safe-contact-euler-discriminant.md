# A same-degree safe insertion can distinguish connectivity beyond the Euler increment

For an R0 occupied prefix, let a vacant next site touch e occupied edge
contacts in c distinct occupied components. Write

```
loop = e-c,
merge = max(c-1,0),
isolate = 1{c=0}.
```

Adding that site changes the graph Betti numbers by
`Delta beta_0=isolate-merge` and `Delta beta_1=loop`. Consequently
`Delta chi=1-e` regardless of the component partition. If the insertion
preserves ambient rank zero, all the new graph cycles are contractible in
the torus; they need not be irrelevant to later essential births.

This gives a direct mechanism discriminator without fitting a surrogate.
Compare two independent uniformly drawn next labels U,V that both preserve
R0 **and have the same e**. They have identical rank and Euler increments.
For e>0, their loop and merger differences satisfy

```
loop(U)-loop(V) = -(merge(U)-merge(V)).
```

For e=0 both differences vanish. A nonzero covariance between the loop
difference and the future birth response therefore isolates attachment
partition information invisible to this immediate rank/Euler summary.
The opposite loop/merger signs in this subset are one exact identity, not
two independent discoveries.

## Readout on the existing next-label forks

For the conditional future response m(Z,u), estimate the same-degree
contrast with

```
1{R0, U and V safe, e(U)=e(V)}
    * (loop(U)-loop(V)) * (Xbar_U-Xbar_V) / 2,
```

where each Xbar averages the two already acquired independent suffixes.
Conditional on the original prefix, its expectation is

```
sum_e pi_safe,e^2 Cov(loop,m | Z, next label safe and contact degree e).
```

There is no division by a sparse empirical class probability. The exact
null `m(Z,u)=f_Z(e(u))` within topology-safe labels makes this contrast zero
for every prefix, even with an arbitrary prefix-specific function f_Z.
Conversely a nonzero population contrast would exclude this particular
rank/Euler-only next-response closure on the sampled prefix population.
A zero contrast would not prove that closure.

The recorded future responses are F1/F2 at the fixed p_ref and their
integrals, equivalently expected K1/K2. Single-orientation contributions
are pooled equally with full original-prefix denominators. This is not an
H4-projected response, a manipulation of loop count alone, a continuum
operator identity, or new independent-prefix evidence. Any error estimate
retains the same twenty original batches per size.

The pre-readout script is `scripts/p334_safe_contact_response.py`, using
the completed e32a8593 conditional tails and the new checkpoint-only
contact coordinates. Numerical results will be reported separately.
