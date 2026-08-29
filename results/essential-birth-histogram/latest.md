# Essential-birth reconstruction from the committed axis L=8 pilot

This is a no-new-production reinterpretation of 100,000 existing samples.

## Archive audit

- every declared SHA-256 checksum matches;
- both marginal totals and the joint total equal `100000`;
- the joint table reconstructs both marginals: `True`;
- all joint integer moments reproduce metadata: `True`;
- support rows with `K_minus>K_plus`: `0`.

## Homology-birth interpretation

The two threshold histograms are the first- and second-essential-birth distributions.
They reconstruct `P0=1-F1`, `P1=F1-F2`, `P2=F2`, and `M=P2-P0` exactly.
At the archived finite root `0.59258424993389151231186158469548399009472832417762`, the recomputed `M` is `3.7371523569099677832098989449490000000000000000000E-50` and the equal-weight birth
mixture CDF is `5.0000000000000000000000000000000000000000000000002E-1`.

The neutral-area identity becomes an exact priority lifetime:

```text
integral P(R=1) dp = E[tau_second-tau_first]
                    = 110129/1300000
                    = 8.47146153846153846E-2 (decimal)
```

The joint archive also determines `E[C]`, `Var(C)`, and `Var(W)` for
`C=(tau_first+tau_second)/2` and `W=tau_second-tau_first`.

## Missing marks

The historical files do not contain the projective winding line `ell`, integral saturation
index, or first/second birth-site local marks. Those quantities require a future stream; they
cannot be reconstructed from marginal or joint endpoint counts.

This analysis reuses the same pilot data and is not new independent evidence.
