# Exact small matching polynomials: axis versus diamond

This note records brute-force exact data from `scripts/exact_matching_polynomial.py`.

The polynomial is the finite Mertens-Ziff matching function written from the wrapping side,

\[
M(p)=\sum_{k=0}^{N} a_k p^k(1-p)^{N-k},
\]

where `a_k` is the integer sum of

`1{black NN wraps} - 1{white NN+NNN wraps}`

over all configurations with exactly `k` occupied sites.

The power-basis coefficients are therefore exact integers. Numerical roots below are diagnostics only.

## Axis torus

### L = 2, N = 4

\[
M(p)=-1+4p^2-2p^4.
\]

Unique physical root:

```text
0.541196100146197...
```

### L = 3, N = 9

\[
M(p)=-1+6p^3-18p^7+18p^8-4p^9.
\]

Unique physical root:

```text
0.586511455112676...
```

### L = 4, N = 16

\[
\begin{aligned}
M(p)={}&-1+8p^4+32p^6-64p^7+172p^8-704p^9\\
&+1104p^{10}-608p^{11}-56p^{12}+128p^{13}\\
&+16p^{14}-32p^{15}+6p^{16}.
\end{aligned}
\]

Unique physical root:

```text
0.590672112331028...
```

The axis roots at these tiny sizes lie below the accepted neighborhood of the infinite-lattice threshold and move upward rapidly.

## Diamond torus

Here the quotient periods are `(L,L)` and `(-L,L)`, so `N=2L^2` and the physical period is `sqrt(2)L`.

### L = 2, N = 8

\[
M(p)=-1+28p^4-48p^5+24p^6-2p^8.
\]

Unique physical root:

```text
0.604563277853507...
```

### L = 3, N = 18

\[
\begin{aligned}
M(p)={}&-1+126p^6-216p^7+180p^8-660p^9+1332p^{10}\\
&-1368p^{11}+1548p^{12}-1980p^{13}+1494p^{14}\\
&-468p^{15}-18p^{16}+36p^{17}-4p^{18}.
\end{aligned}
\]

Unique physical root:

```text
0.594252321168569...
```

The diamond roots at these tiny sizes lie above the infinite-lattice threshold neighborhood and move downward.

## Immediate signal

At the smallest exactly enumerable sizes, **the same matching-root estimator approaches from opposite sides in the two orientations**:

| geometry | L | N | root |
|---|---:|---:|---:|
| axis | 2 | 4 | 0.541196100146197... |
| axis | 3 | 9 | 0.586511455112676... |
| axis | 4 | 16 | 0.590672112331028... |
| diamond | 2 | 8 | 0.604563277853507... |
| diamond | 3 | 18 | 0.594252321168569... |

This is **not yet an asymptotic exponent or spin assignment**. The geometries do not even have matching site counts at the same integer `L`, so a direct average of these rows is meaningless. The relevant comparison must use physical circumference and matched aspect/modular geometry.

Nevertheless, the sign pattern is exactly the qualitative behavior needed for the orientation-projection program: a geometry-dependent leading amplitude is visibly large in the target square-site matching estimator itself.

## Factor / GCD smoke test

SymPy factorization over `Z[p]` leaves each of the five displayed polynomials without a useful common low-degree factor. In particular:

```text
gcd(axis L=2, axis L=3, axis L=4) = 1
gcd(diamond L=2, diamond L=3) = 1
```

This is only a tiny-size negative control. It does **not** show that `p_c` is non-algebraic. It says that the most naive mechanism seen in exactly solvable matching-polynomial examples — one basis-independent low-degree factor already present at every small size — is absent here.

## Next exact targets

1. Reproduce published axis finite matching/critical polynomials as an external regression test.
2. Replace brute-force `2^N` enumeration with a connectivity transfer/dynamic program so diamond `L>=4` becomes accessible.
3. Compute roots on several sheared tori with comparable physical areas/circumferences.
4. Fit the **signed geometry amplitude**, not merely separate threshold sequences.
5. Use modular/physical length consistently before any cancellation claim.
