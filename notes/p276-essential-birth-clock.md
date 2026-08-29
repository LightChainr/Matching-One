# Essential-birth spectroscopy of the threshold-rank archive

## Exact coordinate change

For every Newman--Ziff permutation, let `K1` and `K2` be the first occupation
ranks at which the ambient homology rank becomes at least one and then two.
The existing threshold-rank archive already stores sufficient statistics for

\[
C=\frac{K_1+K_2}{2(N+1)},\qquad
W=\frac{K_2-K_1}{N+1}.
\]

These are not learned latent coordinates:

- `C-1/2` is the complement-odd translation of the two-birth clock;
- `W` is the complement-even lifetime of the unique rank-one interval;
- `E[W]` is the integrated probability of exactly one ambient torus cycle;
- the stored `K1,K2,K1^2,K2^2,K1*K2` sums determine the full `2 x 2`
  covariance of `(C,W)` exactly.

No random field was replayed for this analysis.

## Retrospective result

Ten archived same-`N` orientation pairs were rewritten in `(C,W)`, including
the P43, P49, P50, P57 high-statistics blocks and the P154 variance pilots.
For every pair, both normalized H4 contrasts have a common sign after division
by `Delta cos(4 theta)`:

```text
C / DeltaCos4 < 0 in 10/10 pairs,
W / DeltaCos4 > 0 in 10/10 pairs.
```

The high-statistics P43/P49/P50/P57 subset gives

```text
N^(13/8) DeltaC/DeltaCos4:
  common amplitude = -0.3361312,
  chi2 = 11.0268 / 7.

N^(11/8) DeltaW/DeltaCos4:
  common amplitude = +0.0414431,
  chi2 = 4.4534 / 7.
```

The lifetime power is not identified uniquely by this reuse.  Fixed
`N^-5/4` is almost equally compatible (`chi2=4.4635/7`), whereas `N^-1`
gives `10.9417/7`.  An unweighted log-log diagnostic estimates decay powers
`1.676` for `C` and `1.362` for `W` after excluding the explicitly
underpowered norm-4 variance pilots.

The useful scientific statement is therefore the coordinate split, not a new
claim that `11/8` is an asymptotic exponent.

## Mechanism implication

The global H4 state has at least two topology-defined components:

```text
clock translation C     -- matching/complement odd, clean N^-13/8-like transfer;
rank-one persistence W  -- complement even, slower and distinct finite-size transfer.
```

This supplies a concrete candidate for the secondary even direction seen by
the norm-4 Jordan analysis.  It is more constrained than adding an arbitrary
bulk field: the second coordinate is the exact lifetime of the intermediate
rank-one homology state.

The next stream should record, at `K1`, the primitive rank-one line `ell` and
its Smith index `iota`.  The marked persistence state

\[
(C,W,\ell,\iota)
\]

would then separate continuum clock translation, lifetime/shear, homology
polarization and finite-quotient arithmetic memory without introducing a
target-learned basis.

## Boundary

The coordinate identities and moment reconstruction are exact.  The common
power scores are post-reveal, cross-lineage diagnostics from already correlated
archives; they do not create new independent evidence blocks or establish an
asymptotic exponent.
