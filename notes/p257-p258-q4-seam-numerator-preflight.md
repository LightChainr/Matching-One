# Q=4 seam numerator preflight for Issues 257 and 258

## Outcome

The repository does not yet have an orientation-resolved Q=4 colour-seam H4
kernel.  It does, however, have enough exact machinery to close every semantic
gate around that missing kernel.  A 2-by-2 critical Potts row transfer now
implements the identity and `(01)` seams, exact central projectors for the
`[4]` singlet and `[2,2]` sectors, two independent nonzero equivariant transfer
insertions, and an independent direct spin-configuration enumeration.

The unnormalized sector numerators reproduce the frozen #257 character rule:

```text
                         identity      (01) seam       ratio
singlet [4], T^2           15904          15904           1
[2,2], T^2                    32              0           0
singlet [4], d_log(v)T^2    63616          63616           1
[2,2], d_log(v)T^2            128              0           0
```

The two witnesses are deliberately different, so the result is not a replay
of the pair-space character calculation.  They are full row-transfer traces
`Tr(P_lambda U_g K)` on a 16-state row space.  Central-projector commutation
with the transfer makes Schur factorization executable:

```text
Tr(P_lambda U_g K) / Tr(P_lambda K)
    = chi_lambda(g) / dim(lambda).
```

This is the positive control needed before introducing a physical H4 kernel.

## The #258 normalization cross-check

The full tiny-torus partition sums are seam dependent:

```text
Z_identity = 33024,
Z_(01)     = 21568,
Z_(01)/Z_identity = 337/516.
```

The row-transfer trace and an independent enumeration of all `4^4=256` spin
configurations give the same integers.  Consequently, a normalized singlet
expectation would have seam ratio

```text
(N_(01)/Z_(01)) / (N_identity/Z_identity) = 516/337,
```

not one.  Restoring the partition factor gives

```text
[E_(01)/E_identity] [Z_(01)/Z_identity] = 1.
```

This is the finite-volume bridge between #257 and #258: representation
characters control the observable numerator, while the measure/partition
normalization is a separate object.  They must not be conflated.

## Exact boundary

The exact result applies to any S4-equivariant insertion inside the projected
isotypic transfer sector.  The `T^2` and `d_log(v)T^2` rows prove the seam and
normalization plumbing with nonzero numerators.

The log-v insertion is not a pure continuum H4 operator.  On a fixed square
bond stencil, scalar and fourth-harmonic microscopic components alias unless
an orientation-resolved kernel is supplied.  Therefore this commit does not
claim a V22 measurement and does not start production.

## Frozen minimal production schema

The next runner must add only one missing component: an orientation-resolved
H4 transfer insertion with declared handedness.  For one fixed torus cycle it
must emit raw unnormalized values

```text
Z_identity, Z_transposition,
N_H4_singlet_identity, N_H4_singlet_transposition,
N_H4_[2,2]_identity, N_H4_[2,2]_transposition.
```

The primary scores are fixed at `1` and `0`; no exponent is fitted.  If the
implementation emits normalized `E_g=N_g/Z_g`, the scorer must reconstruct
`N_g` or equivalently multiply the expectation ratio by
`Z_transposition/Z_identity`.  Production remains blocked until:

1. the tiny partition sums in this artifact reproduce exactly;
2. both projectors commute with the zero-source transfer;
3. both identity-seam H4 numerators are nonzero;
4. the H4 kernel's orientation and handedness are explicit and it is not the
   log-v witness used here.

No large runner or general FK framework is needed.

## Scientific layers

**Exact.** The tiny Potts transfer, direct spin sum, sector traces, `1/0`
ratios, and normalization restoration are rational/integer identities.

**Mechanism inference.** One colour transposition is enough to distinguish a
singlet spin-4 numerator from a `[2,2]`-charged numerator without a scaling fit.

**Exploratory conjecture.** If a future controlled `[2,2]` H4 numerator is
nonzero, it supplies the V22-like channel required to compare the local
four-arm mark with the globally singlet thermal-Q4 channel.  This preflight
does not assert that overlap.

## Reproduction

```bash
python3 scripts/exact_q4_seam_numerator_preflight.py \
  --output results/q4-seam-numerator-preflight/latest.json
python3 -m unittest discover -s tests \
  -p 'test_exact_q4_seam_numerator_preflight.py'
```

