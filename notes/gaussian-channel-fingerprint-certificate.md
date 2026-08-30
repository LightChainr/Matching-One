# Exact four-channel Gaussian fingerprint certificate

This certificate addresses the bounded algebraic part of parent issue #64.
It turns the four normalized-P4 multiplier targets into rational fingerprints,
without introducing decimal tolerances.

If a channel has radial exponent `k/8`, its normalized multiplier is
`Q^(-k/8)`.  Raising to the eighth power gives the exact rational fingerprint
`Q^(-k)`.  For the frozen exponents

| channel | `k` |
| --- | ---: |
| `P4_S` | 8 |
| `P4_Dprime` | 5 |
| `P4_D` | 13 |
| `P4_Sprime` | 10 |

all four fingerprints are distinct at both `Q=2` and `Q=5`.  The derivative
pairs also obey exact relations

`fingerprint(P4_Dprime) / fingerprint(P4_S) = Q^3`

and

`fingerprint(P4_Sprime) / fingerprint(P4_D) = Q^3`.

For the norm-5 same-radial H12 adversary, the normalized angular factor is
`-1679/625`.  The certificate records the negative sign separately and freezes
the exact eighth-power magnitude for `P4_D` and `P4_Sprime`.

## Boundary

This closes only the exact non-aliasing and sign/amplitude arithmetic.  It does
not read production or covariance data, score residuals, or choose a physical
channel.  Parent issue #64 therefore remains open.
