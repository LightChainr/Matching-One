# Exact parity and harmonic algebra for composite spin-4 corrections

Status: independent selection-rule check for Issue 58. This does not score production data.

## Declared generators

The oracle uses exactly the generator set proposed in the issue:

| generator | matching parity | spin | length exponent | status |
|---|---:|---:|---:|---|
| `T4` | -1 | 4 | `13/4` | observed thermal-family candidate |
| `I4` | +1 | 4 | `2` | observed identity-family candidate |
| `V4` | +1 | 4 | `8/3` | conditional parity assignment |
| `S0` | +1 | 0 | `2` | optional even scalar |

Every spin-4 factor contributes `cos(4 theta)`. The leading matching-odd H4 term is one `T4`,
so a monomial's relative exponent is its total length exponent minus `13/4`; the accelerated
axis-root exponent is `w=4+q`.

## Exact angular calculation

Set `z=exp(4 i theta)`. Then

```text
cos(4 theta) = (z+z^-1)/2.
```

The script multiplies this two-term Laurent polynomial with rational coefficients. Combining
the `z^k` and `z^-k` terms gives the real cosine harmonic `H_(4k)`. This proves two counting rules:

1. matching oddness requires an odd number of `T4` factors;
2. H4 support requires an odd total number of spin-4 factors.

Consequently `T4*I4` is matching odd but contains only H0 and H8. A single even spin-4 insertion
cannot correct the H4 channel.

For three spin-4 factors the exact expansion is

```text
cos(4 theta)^3 = (3/4) cos(4 theta) + (1/4) cos(12 theta).
```

Thus each named cubic product has both H4 and H12 support at the elementary angular level:

| monomial | q | w | H4 | H12 |
|---|---:|---:|---:|---:|
| `T4*I4^2` | `4` | `8` | `3/4` | `1/4` |
| `T4*I4*V4` | `14/3` | `26/3` | `3/4` | `1/4` |
| `T4*V4^2` | `16/3` | `28/3` | `3/4` | `1/4` |

The `1/3` H12/H4 ratio is exact for one aligned monomial before response-tensor or continuum mixing.
The robust output for the repository is joint harmonic support, not an unconditional measured ratio.

## All-degree q=3 no-go inside the declared semiring

Let `(t,i,v,s)` be the nonnegative multiplicities of `(T4,I4,V4,S0)`. Matching oddness requires
positive odd `t`. Relative to one leading `T4`,

```text
q = 13(t-1)/4 + 2(i+s) + 8v/3.
```

Multiplying by 12 gives

```text
12q = 39(t-1) + 24(i+s) + 32v.
```

For `q=3`, the right side must equal 36. If `t>=3`, the first term alone is at least 78. If
`t=1`, division by four gives

```text
6(i+s) + 8v = 9,
```

which is impossible because the left side is even. Therefore no analytic monomial at any degree
in this declared generator semiring produces `q=3` (`w=7`). This is stronger than a finite table,
but it remains conditional on the generator set being complete.

## q=6 is not exponent-unique

The ordinary next thermal spin-4 quasiprimary has `q=6`, `w=10`. The same exponent is also generated,
when the optional scalar exists, by analytic composites including

```text
T4*S0^3,
T4*I4^2*S0.
```

Therefore `w=10` alone does not uniquely identify the ordinary thermal tower. Harmonic support,
conditional-scalar controls, and independent amplitudes are needed. The generated artifact lists all
matching-odd H4 collisions through frozen total degree five.

## Evidence boundary

Exact statements:

- matching parity multiplication;
- Laurent/Fourier coefficients and H4/H12 support;
- rational exponent sums;
- the all-degree `q=3` exclusion inside the declared generator semiring.

Conditional or unproved statements:

- the even parity of `V4`;
- existence and nonzero coupling of `S0`;
- completeness of the four generators;
- nonzero amplitudes of any composite;
- survival of the elementary `1/3` ratio after continuum mixing.

The parent issue remains open for held-out annihilator scores, Gaussian residual power, and a measured
subleading H12 channel.
