# Rank-plane crosswalk: clocks, projectors and typed insertions

## Exact basis

```text
clock odd/even:       c=C-1/2                 W=(K2-K1)/(N+1)
state odd/even:       A=P2-P0                 E=P0+P2=1-P1
birth odd/even:       D=f12-f01=E'            S=f01+f12=A'
historical notation:  D_hist=A/2, S_hist=E/2, S_hist'=D/2, D_hist'=S/2
endpoint mixtures:    F2=(A+E)/2              f12=(S+D)/2
```

Complement/reversal makes `c,A,D` odd and `W,E,S` even. Therefore
the historical `K2`/`F2` cancellation cannot be the even `W/S`
projector: `F2` is exactly a mixed-parity endpoint.

## Archive result

The joined high-statistics streams cover `N=[130, 145, 170, 185, 265, 290, 325, 425]`.
The median odd shares are `0.808`
in the clock plane, `0.607` in the
state plane, and `0.819` in the
birth-density plane. The median K2 cancellation fraction is
`0.258` at the state layer,
versus `0.617` for its clock
endpoint and `0.639` for its
density endpoint. The strong cancellation is therefore localized to
the state value at the intrinsic center; it is not a universal K2 null.

Each dataset in the JSON carries a complete `14 x 14` same-batch
covariance in this basis. The fixed-center matrix is an exact covariance
of the archived batch estimators. The intrinsic-center matrix adds the
displayed first-order M-estimator influence of the fitted center.

## What is and is not reconstructible

The production moments recover event-level `Cov(C,W)`, while marginal
histograms recover all endpoint CDF/density means. Common batch IDs
recover their estimator covariance. Event-level nonlinear cross moments
need the joint `K1,K2` histogram; `ell`, `iota`, `chi4(ell)D`, and `qJ`
cannot be recovered from the old production files.

## Next acquisition

Among committed single-size geometries, `P50-N145`
(`N=145`, {'first': [12, 1], 'second': [9, 8]}) maximizes the
outcome-free `|Delta cos4|/sqrt(N)` mark-acquisition proxy. The more
informative minimal radial campaign is the frozen q2 edge
`N65 -> N130`: its `(1-i)` multiplier flips `chi4(ell)` exactly and is
the cheapest way to separate clock translation, rank-one lifetime and
line polarization across scale.

Store a sparse joint histogram of
`(K1,K2,ell_u,ell_v,iota,count)` per aligned batch/orientation. Add
`q`, `J_D4`, and `q*J_D4` batch sums only when the connected #275
coupling is desired.

## Scientific card

1. MECHANISM SPACE: K2 cancellation is resolved as A_top/E_top cross-parity interference; the exact eigenplanes are (C,W), (A,E), and (D_birth,S_birth).
2. NOT PROVED: archives do not establish a one-dimensional Q4 field, an ell/iota continuum overlap, or an asymptotic exponent for W/S_birth.
3. OBSERVER-SECTOR-SOURCE-GEOMETRY: A_top/E_top | odd/even Alexander sectors | chi4(ell)D/S sources | Gaussian same-N orientation pairs.
4. DEPENDENCY GROUP: all reconstructed contrasts and covariances are views of the same threshold-rank streams, not independent evidence rows.
5. UPWEIGHT OBSERVATION: acquire sparse joint (K1,K2,ell,iota) on the q2 N65->N130 sign-flip edge, with qJ_D4 cross moments when the connected #275 coupling is targeted.

## Claim boundary

The basis and parity statements are exact. Archive covariance is
same-stream reuse; the intrinsic-center extension is first-order.
Sign coherence and odd-leading language are mechanism diagnostics,
not independent evidence or an asymptotic field identification.
