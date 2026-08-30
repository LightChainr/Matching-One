# P250 projective-leg operator replacement

Status: the theory-selected projective-leg row passes the frozen 2k
propagation/interface gate.  This is the first P250 insertion in the current
sequence to resolve every charged pair denominator at two separations.

## Why this is a different operator

The failed R1--R4 family asked a local pivotal event to carry information from
an annular landing stencil.  The new row attaches the root directly to a
rank-one wrapping component:

```text
L(x) = 1{x in a black NN rank-one component}
     - 1{x in a white matching rank-one component}.
```

It is therefore mesoscopic by construction.  If two roots lie in the same
rank-one component they share an exact contribution.  The Z5 charge still
comes from the same five-fiber DFT, so `G_r(d)`, `C113` and `C122` retain the
existing cover-character interface.

The scalar leg row was chosen before sampling.  A primitive-line `chi4`
version is exact and remains available as a secondary typed row, but was not
allowed to win or lose through angular cancellation.

## Exact smoke

For both norm-5 hands, a Manhattan representative of the first primitive
period gives a rank-one wrapping component.  At four separated roots its
black marks are exactly

```text
[+1,+1,+1,+1],
```

and the complementary white-matching construction gives

```text
[-1,-1,-1,-1].
```

The component line has gcd one.  Its physical `chi4` is invariant under sign
and an exact unimodular period-basis shear.  The maximum numerical residual in
the declared Z5 deck-character shift is `2.49e-16`.

## Frozen 2k result

Seed `25033433720260830`, counters `[0,2000)`, 40 batches:

| d | weakest of four pair z scores | pair gate | cubic covariance | descriptive cubic support |
|---:|---:|---|---|---:|
| 1 | 5.527 | pass | nondegenerate | `15.204/8`, p `0.0553` |
| 2 | 2.900 | pass | nondegenerate | `28.170/8`, p `0.000443` |
| 3 | 0.709 | fail | nondegenerate | `25.825/8`, p `0.00112` |

The promotion rule uses only the first three columns.  It passes at `d=1,2`.
No phase or closure statistic was computed.  The cubic support values show
that the joint interface is not merely numerically invertible, but they did
not select the operator.

The pair means also display the intended propagation scale: they are roughly
`0.0066--0.0084` at `d=1` and `0.0030--0.0045` at `d=2`, before becoming
under-resolved at `d=3`.  This is qualitatively different from the local-H4
radius scan, where no radius resolved all four denominators beyond contact.

## Decision

Advance this operator, not the smoke data.  Freeze a fresh counter block using
the identical projective-leg definition and score `d=1,2` support before any
phase model.  Keep `d=3` as an optional tail diagnostic, not a required gate.

The passing row is nonlocal/topological, so it cannot silently replace a local
primary-field claim.  It establishes a usable charged projective-sector
observable; continuum identity and OPE normalization remain separate tasks.
