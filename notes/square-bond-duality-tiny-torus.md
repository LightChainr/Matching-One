# Square-bond self-dual control: exact tiny-torus conventions

Status: C5 finite identity / C1 method freeze for issue #42. Not an orientation-amplitude result.

## Why this slice exists

Issue #42 asks whether the two-sector picture

```text
matching-even  S4 ~ L^{-2}
matching-odd   D4 ~ L^{-13/4}
```

survives square-bond percolation at the exact self-dual point `p=1/2`. Production Gaussian runs are still required for the orientation projector. They cannot start until primal/dual wrapping channels and the even/odd combinations are fixed by identities that do not depend on a Monte Carlo seed.

## Dual transport, not bit-complement

On an `L x L` square torus each of the `2L^2` primal bonds is paired with the unique dual bond that crosses it. Under the vertex identification already used by `square_bond_pairs`, every dual displaced edge equals some primal displaced edge.

Geometric dual transport is the bijection `T` that occupies, for each vacant primal bond, the identified primal edge equal to that bond's dual edge. Naive bit-complement of the same mask occupies the vacant primal edges instead. Exhaustive enumeration shows that bit-complement does **not** swap primal and dual wrapping (138 failures on `L=2`, 147560 on `L=3`). Treating complement as the duality involution is a convention error.

`T` itself is a bijection. For every enumerated configuration and every registered channel (`cross`, `both`, `either`, `direction_0`, `direction_1`) it swaps primal and dual wrapping, so

```text
S = (R_primal + R_dual) / 2
D =  R_primal - R_dual
```

satisfy `S(T(ω))=S(ω)` and `D(T(ω))=-D(ω)`. At `p=1/2` every configuration is equiprobable, therefore

```text
E[D(p=1/2)] = 0
```

for every finite `L` and every channel. Equivalently, primal and dual wrapping are equal in law: the dual graph is the same square grid and dual occupation is i.i.d. Bernoulli(`1/2`). This is exact, not asymptotic. It is the identity-level content of the statement "the duality-odd central amplitude vanishes at the self-dual point".

The even combination `S` is not forced to zero. Locked rationals:

```text
L=2, 256 configurations
  cross        E[S] = 69/256
  both         E[S] = 73/256
  either       E[S] = 187/256
  direction_0  E[S] = 65/128
  direction_1  E[S] = 65/128

L=3, 262144 configurations
  cross        E[S] = 18865/65536
  both         E[S] = 20785/65536
  either       E[S] = 46671/65536
  direction_0  E[S] = 527/1024
  direction_1  E[S] = 527/1024
```

The `either = direction_0 + direction_1 - both` channel algebra holds exactly on both tables.

## What this does not establish

- orientation contrast `P4[S]` or `P4[D']` on same-`N` Gaussian pairs;
- the candidate exponents `L^{-2}` and `L^{-13/4}`;
- a continuum duality-parity assignment for thermal Q4;
- any statement about square-site `p_c`.

Same-`N` orientation pairs in the current Gaussian catalog begin at `N=65`. Tiny exact quotients do not supply two D4-inequivalent orientations of equal area, so they cannot score the projector.

## Remaining production content of #42

1. Freeze sample counts from a pilot on `N=65` `(8,1)` vs `(7,4)` at exact `p=1/2`.
2. Measure `P4[S]` (duality-even anisotropy) and `P4[D']` (odd-sector slope) under common random numbers, using the primal/dual/`T` conventions frozen here.
3. Compare against the frozen square-site/matching amplitudes only after the observable analogue and torus modulus are matched.

The exact identity above already says that the *center value* of `D` cannot carry an orientation signal at `p=1/2`. Any surviving odd-sector signal must live in derivatives or in a full-curve reconstruction around `1/2`.

## Oracle

```bash
python3 scripts/square_bond_duality_exact.py --L 2
python3 scripts/square_bond_duality_exact.py --L 3 --json
```

Regressions live in `tests/test_square_bond_duality_exact.py`.
