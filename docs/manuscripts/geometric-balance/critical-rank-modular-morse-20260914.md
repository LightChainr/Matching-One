# Critical torus rank-shape geometry in modulus space

Date: 2026-09-14

Status: exact Q=1 formula + independent high-precision controls + a sharply reduced theorem target for #781. No global minimum theorem is claimed here.

## 1. Critical rank law

For critical Q=1 FK percolation on a torus of modulus `tau`, Arguin/Pinson homology sums give

```text
P0(tau)=P2(tau).
```

Use the unit-area Gaussian lattice sum

```text
Theta_g(tau)
 = sqrt(g) / (sqrt(Im tau) |eta(tau)|^2)
   * sum_{m,n in Z} exp[-pi g |m tau-n|^2 / Im tau].
```

The aggregate trivial/cross probability is

```text
P0(tau)
 = 1/2 [Theta_{8/3}(tau)-Theta_{2/3}(tau)].              (1.1)
```

The intrinsic even rank-shape coordinate is

```text
c_*(tau) = P1/(2 P0) = 1/(2 P0)-1.                      (1.2)
```

Thus minimizing `c_*` is equivalent to maximizing `P0`.

## 2. A useful theta-ratio reduction

Let

```text
theta(alpha;tau)
 = sum_{m,n} exp[-pi alpha |m tau+n|^2/Im tau],
D(tau)=sqrt(Im tau)|eta(tau)|^2.
```

Two standard identities are

```text
theta(alpha;tau)=alpha^-1 theta(alpha^-1;tau),
D(tau)=sqrt(6)/4 [2 theta(6;tau)-theta(3/2;tau)].
```

Substituting into (1.1) gives the exact scalar ratio

```text
P0(tau)
 = [(4/3) theta(8/3;tau)-theta(3/2;tau)]
   / [2 theta(6;tau)-theta(3/2;tau)].                    (2.1)
```

This is a more precise proof target than a generic statement that "hexagonal lattices extremize theta functions". Existing results of Montgomery and Luo--Wei concern individual theta functions, selected differences, and ratios; (2.1) is a signed affine ratio and is not silently covered without checking theorem hypotheses.

Relevant boundaries:

- Luo--Wei, arXiv:2203.00264, proves hexagonal minimization for `theta(alpha)-beta theta(2 alpha)` in its stated parameter range.
- Luo--Wei, arXiv:2605.07580, classifies broad theta/Epstein ratio extrema and includes the Dedekind-eta/free-boson determinant as an application.

Neither citation is used here as a proof of the #781 inequality.

## 3. Independent Morse controls

At the equianharmonic modulus

```text
tau_hex = 1/2 + i sqrt(3)/2,
```

80-digit Gaussian summation gives

```text
P0 = 0.316053412701401686065271069685098217779...
c  = 0.582011077578161892499173198398907452715...
```

The numerical gradient vanishes to the working precision. In ordinary coordinates `(Re tau, Im tau)`, the Hessian is

```text
H_hex = 1.01315266106979019658579486489... * I.          (3.1)
```

Hence the hexagonal point is a strict nondegenerate local minimum of `c_*` (equivalently a strict local maximum of `P0`). The scalar Hessian is consistent with the order-three modular stabilizer.

At the square modulus `tau=i`,

```text
P0 = 0.309526275429831322769693633352898776545...
c  = 0.615371746084117174271794028861523956059...
```

and

```text
H_square = diag(
 -0.482066947400697785573710130065...,
  1.92866506353063892452099278454...
).                                                        (3.2)
```

So the square point is a saddle, not a competing local minimum.

These are numerical certificates, not a global theorem.

## 4. Source-zero collision locus

At critical balance the bounded rank-source partition function is

```text
Z_X(s;tau) = [c_*(tau)+cosh s]/[c_*(tau)+1].             (4.1)
```

Therefore:

- `c<1`: nearest zeros are purely imaginary, `s=+- i acos(-c)`;
- `c=1`: the two zeros coalesce at `s=i pi` into a double zero;
- `c>1`: the pair leaves the imaginary axis, `s=i pi +- arcosh(c)`.

Thus

```text
c_*(tau)=1                                                 (4.2)
```

is an exact modular source-zero collision locus. It is a zero-geometry bifurcation of a bounded topological-source polynomial, not a thermodynamic phase transition and not a Jordan diagnostic.

For rectangular `tau=i r`, the previously computed crossing is

```text
r_* = 1.78782935267965703762908251131828066576...
```

(and its modular reciprocal).

A direct scan over `0<=Re tau<=1/2` shows the upper collision branch is extremely close to a horocycle. On eleven equally spaced `x` values, a first harmonic

```text
y_*(x) ~= y0 + A [1-cos(2 pi x)]
```

has

```text
y0 = 1.7878293526796570376290825...,
A  = 1.81122997e-4,
max grid residual ~= 1.46e-7.
```

This is only a numerical description.

## 5. Reduced #781 theorem target

A proof of the conjectured global minimum can be organized around the exact statement

```text
c_*(tau) >= c_*(tau_hex)
```

or equivalently the signed theta inequality (2.1) with maximum at `tau_hex`.

A practical proof architecture is:

1. use modular invariance to classify the square and hex stationary points;
2. prove the hex Hessian scalar is positive analytically;
3. prove monotonicity on the fundamental-domain boundary arcs;
4. exclude additional interior critical points using theta heat-equation / lattice-moment identities;
5. only then invoke the closest applicable Luo--Wei theta-ratio inequalities.

The current numerical evidence strongly supports the hexagonal conjecture but does not replace steps 3--4.

## 6. Claim boundary

- Exact: (1.1), (1.2), (2.1), (4.1)--(4.2).
- Certified numerical controls: (3.1)--(3.2), rectangular collision and coarse collision curve.
- Conjecture: global hex minimum and the simple topology of the `c<1` core in the modular fundamental domain.
