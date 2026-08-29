# PR267 Target 1: leave the q-only algebra with one Euler coordinate

## Decision

Do not start the old N325/N425 stream. Its primary connected row would still
be a function of the three-state ambient rank variable and a root-independent
S/D gate, hence contact-closed by the finite algebra in `b31581a`. The bare
primitive-line scalar transfer and the frozen local-D4 UV selector have also
failed their held-out tests (`507571d`, `eea85c4`).

The minimal external observer is instead

```text
O_ext(A) = C_black^NN(A) - C_white^matching(A^c) - q(A).
```

The configuration-level Euler identity already implemented on the integer
period tori gives the second exact form

```text
O_ext(A) = V(A) - E_NN(A) + F0(A).
```

Thus this one number has two useful readings. The component form is almost
free inside the two union-find paths. The cell form identifies the new
coordinate as a bulk Euler/thermal scalar. At fixed microcanonical `k`, its
fluctuation is exactly `-(E-E_mc)+(F0-F0_mc)`.

## Why this escapes the no-go

The old topology-only observer algebra is `sigma(q)` with
`q in {-1,0,+1}`. The index-9 control contains an explicit fixed-`k` pair
with equal `q` and unequal `O_ext`; axis L2 and Gaussian index 5 separately
certify the component/cell and path identities.
Therefore no polynomial or arbitrary function of `q` represents `O_ext`.
The all-order contact theorem cannot reduce `O_ext*J` to q-only gate moments.
It may still vanish dynamically, but that is now an empirical field question,
not an algebraic tautology.

Under primal/matching complement, both `C_black-C_white` and `q` reverse, so
`O_ext` is matching-odd. Consequently `O_ext*J_D4` is matching-even and is the
primary allowed bridge; `O_ext*J_S4` is an odd companion/null channel. The
old `q*J_D4` stays in every row only as an exact contact control.

## One stream, several decisions

For every pre-insertion size and aligned batch the runner now retains

```text
O_ext, O_ext^2,
O_ext J_D4, O_ext J_S4,
J_D4 conjugate(J_S4), |J_S4|^2.
```

The last two are same-next-site Horvitz-score Gram entries. They are not
mislabelled as products of two independently summed site sources. Together
they permit covariance, partial-covariance, and two-source Gram projections
without another simulation. `K1/K2`, typed births, `ell/iota`, line/landing
H4, S/D, and qJ remain in the same stream.

## Candidate comparison

1. **Bulk Euler/Betti residue — selected.** It is exact, scalar, outside q,
   matching-typed, and costs two maintained component counters rather than a
   scan of the configuration.
2. **Macroscopically separated local-H4/arm insertion — reserve as the second
   coordinate.** It is field-specific and genuinely two-point, but requires a
   frozen far-anchor orbit on every integer-period quotient. Reusing the old
   one-point radius convention would simply repeat the rejected selector.
3. **Charged seam/winding character — keep as a charged-sector experiment.**
   It can escape the neutral algebra, but only after fixing a seam gauge,
   rotation action, and per-component charge update.

The selection is not a claim that the bulk observer must couple. It is a
claim that it is the cheapest exact coordinate for which the coupling is not
already predetermined by the rank-gate algebra.

## Direction-character authorization gate

Main commit `83e98fc` adds a second exact no-go: a direction-only readout on
one C4 orbit aliases scalar and spin four because `exp(-4 i theta)` is
constant on that orbit. This design passes that gate in a specific way:

- `O_ext` is deliberately a scalar Euler coordinate, not a direction-only
  observable relabelled as spin four;
- the H4 source retains typed internal complex information
  `chi4(P ell)`, the primitive line, and separate axis/diagonal landing fields;
- the measured bridge is scalar `O_ext` times complex `J_D4/J_S4`.

The reserved distant local-H4 coordinate is **not authorized** as a future
single-orbit scalar. It must retain separate axis and diagonal C4 orbits, whose
response matrix `[[1,1],[1,-1]]` has determinant `-2`, or carry equivalent
typed/internal complex edge data. This is a pre-production gate, not a
post-hoc diagnostic.

## Production boundary

The code and schema are production-compatible, but the N325 `(17,6)/(18,1)`
and N425 `(16,13)/(19,8)` jobs remain forbidden in this commit. Production
starts only after the external-observer scorer and geometry-pair contract are
frozen. The tiny preflight uses the existing index-9 topology backend and
checks all microcanonical path products coefficient by coefficient.

## Local contact split and revealed pilot

The first N65/N130 100k run found a very strong nonzero `O_ext J_D4`, but its
complex transfer phase was indistinguishable from the `O_ext J_S4` control.
That is exactly the pattern a root-local contact can mimic. Before authorizing
production, the stream was therefore extended by one frozen nuisance and no
new field family:

```text
h_x(A) = 1_A(x) - 1_A(x)1_A(x+e1) - 1_A(x)1_A(x+e2)
         + 1_{A empty on the unit face anchored at x},
O_near = sum_{max(|dx|,|dy|)<=2} h_{root+(dx,dy)},
O_far  = O_ext - O_near.
```

The window is a translation-invariant D4 scalar containing both axis and
diagonal anchors. The stream stores `O_near`, its square, `O_ext O_near`, and
both complex source products, so `O_far` and all relevant covariance entries
are reconstructed exactly without another pass.

With independent seeds and 100 batches at each size, the revealed pilot gave:

| channel | N65 | N130 | amplitude transfer | phase transfer |
|---|---:|---:|---:|---:|
| connected `O_ext J_D4` | `-6.720+1.937i` | `-11.238-3.300i` | 1.675 | 0.566 |
| connected `O_near J_D4` | `-1.807-0.025i` | `-0.330-0.380i` | 0.278 | 0.842 |
| connected `O_far J_D4` | `-4.913+1.962i` | `-10.908-2.920i` | 2.135 | 0.642 |
| connected `O_ext J_S4` | `1.289-0.377i` | `1.612+0.464i` | 1.249 | 0.565 |

All complement counters and exact Russo residuals were zero. `O_far J_D4`
remains overwhelmingly nonzero at both sizes, while the near term shrinks and
does not account for the total. Thus the Euler observer survives the contact
gate. The pilot selects only the following minimal production model:

1. primary: the two complex `O_ext J_D4` vectors;
2. mechanism gate: the two complex `O_far J_D4` vectors are nonzero;
3. nuisance: `O_near J_D4` with the within-size covariance block;
4. frozen contact null: the N325-to-N425 transfer phase of external JD equals
   that of external JS, with amplitudes left free;
5. transfer amplitudes/phases are reported, but two sizes do not define an
   asymptotic exponent.

The frozen scorer is `scripts/score_external_observer_transfer.py`; it uses
the full within-size delete-one covariance and treats the two production
counter/seed groups as independent.

## Scientific card

1. **Mechanism-space change:** Target 1 moves from q-contact response to a bulk Euler scalar times the rank-birth H4 source.
2. **Not proved:** no Q4-epsilon identity or exponent is assigned before a held-out score.
3. **Observer / sector / source / geometry:** Euler residue / matching-odd scalar / J_D4 and J_S4 / integer-period paired tori.
4. **Dependency group:** one common permutation path supplies state, source, cross-products, and both Gram roots.
5. **Upweighted observation:** freeze and score `Cov(O_ext,J_D4)` jointly with the O-JS null and Gram-conditioned source plane.
