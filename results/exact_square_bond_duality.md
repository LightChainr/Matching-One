# Exact square-bond duality identities on tiny tori

Source: `scripts/square_bond_duality_exact.py`.
Claim level: C5 finite identity. Not an orientation-amplitude archive.

At `p=1/2` geometric dual transport `T` swaps primal and dual wrapping, so
`E[D]=E[R_primal-R_dual]=0` in every channel. Naive bit-complement is not
that map. The even combination `S=(R_primal+R_dual)/2` is not forced to zero.

## L = 2, 8 bonds, 256 configurations

| channel | E[S] | E[D] |
|---|---|---|
| cross | 69/256 | 0 |
| both | 73/256 | 0 |
| either | 187/256 | 0 |
| direction_0 | 65/128 | 0 |
| direction_1 | 65/128 | 0 |

```text
naive_complement_swap_failures = 138
geometric_dual_transport_failures = 0
T_bijective = true
passed = true
```

## L = 3, 18 bonds, 262144 configurations

| channel | E[S] | E[D] |
|---|---|---|
| cross | 18865/65536 | 0 |
| both | 20785/65536 | 0 |
| either | 46671/65536 | 0 |
| direction_0 | 527/1024 | 0 |
| direction_1 | 527/1024 | 0 |

```text
naive_complement_swap_failures = 147560
geometric_dual_transport_failures = 0
T_bijective = true
passed = true
```

On both sizes, `E[either] = E[direction_0] + E[direction_1] - E[both]`.

## Boundary

Same-N Gaussian orientation amplitudes for issue #42 are not inferred from these tiny exact tori. The oracle fixes the primal/dual convention so exploratory or confirmatory larger runs can use a common map.
