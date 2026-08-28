# Minimality of the N=1105 four-orientation Gaussian torus

The four-angle `H0/H4/H8/H12` design at `N=1105` is not just a convenient example. It is the **smallest primitive Gaussian square torus that has four D4-inequivalent microscopic orientations at one fixed site count**.

## Lemma

Let

\[
N=a^2+b^2,
\qquad \gcd(a,b)=1,
\]

and identify Gaussian representatives under square-lattice `D4` symmetry (units and complex conjugation). If `k` is the number of distinct prime divisors of `N` congruent to `1 mod 4`, then the number of primitive D4 orientation orbits is

\[
2^{k-1}
\]

for `k>=1`.

Consequently, four distinct primitive orientations require `k>=3`. The smallest possible site count is therefore

\[
\boxed{N_{\min}=5\cdot13\cdot17=1105}.
\]

Its four first-octant representatives are exactly

```text
(33,4), (32,9), (31,12), (24,23).
```

## Proof sketch

Suppose `a^2+b^2=N` is primitive.

1. A prime `p=3 mod 4` cannot divide `N`. If it did, `a^2=-b^2 mod p`; since `-1` is not a quadratic residue mod `p`, this forces `p|a` and `p|b`, contradicting primitivity.
2. `4` cannot divide `N`: squares mod 4 are `0,1`, and a sum divisible by 4 forces both `a,b` even.
3. Hence

   \[
   N=2^\epsilon\prod_i p_i^{e_i},
   \qquad \epsilon\in\{0,1\},\quad p_i=1\pmod4.
   \]

4. In `Z[i]`, each `p_i` splits as `pi_i * conjugate(pi_i)`. For a primitive Gaussian integer `a+ib`, all `e_i` copies must be assigned wholly to one member of that conjugate pair. If both factors divide `a+ib`, then the rational prime `p_i` divides both real and imaginary parts.
5. Thus every distinct `p_i=1 mod 4` supplies one binary conjugation choice. After quotienting the global conjugation and units (`D4`), the number of orientation orbits is `2^(k-1)`.

To obtain at least four orbits, `k>=3`. Exponents do not increase `k`, so the minimum uses the three smallest distinct `1 mod 4` primes: `5,13,17`.

The executable finite verification is `scripts/minimal_four_orientation_gaussian_torus.py`; it checks the formula against brute-force primitive representations for every `1<=N<=1105` and verifies the four `N=1105` representatives.

## Research consequence

This makes the N=1105 experiment an information-theoretically natural endpoint of the exact same-N angular program:

- two-orientation sizes can separate one angular direction from a scalar only after truncation assumptions;
- `N=1105` is the first primitive square torus where four simultaneous orientations exist;
- four orientations are exactly enough to invert the finite harmonic basis `H0,H4,H8,H12` used in PR #77.

It does **not** imply that H16 and higher harmonics vanish. The existing four-angle projector remains a finite harmonic decomposition through H12, not an all-orders scalar projector.
