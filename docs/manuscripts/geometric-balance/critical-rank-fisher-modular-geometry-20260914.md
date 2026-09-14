# Fisher geometry of the critical self-dual torus rank law

Date: 2026-09-14

Status: exact information geometry for the Q=1 critical rank family. This note connects the rank-simplex Fisher sphere on #773 to the modular rank-shape programme #781.

## 1. Critical family is a one-dimensional Fisher meridian

At critical Q=1 on any torus modulus,

```text
P0=P2=a(tau),
P1=1-2a(tau).
```

The multinomial Fisher metric on the simplex is

```text
ds^2 = sum_j dP_j^2/P_j.
```

Along the self-dual critical line this reduces exactly to

```text
ds^2
 = [2/a + 4/(1-2a)] da^2
 = 2 da^2/[a(1-2a)].                                    (1.1)
```

Put

```text
2a = sin^2 phi,     0<=phi<pi/2.
```

Then

```text
ds = 2 dphi.                                               (1.2)
```

The rank-one cusp `(P0,P1,P2)=(0,1,0)` is `phi=0`. Therefore the Fisher geodesic distance from that cusp is

```text
D_F(tau)
 = 2 asin sqrt(2 P0(tau)).                                (1.3)
```

This is the ordinary spherical distance under the square-root embedding of the probability simplex.

## 2. Exact relation to the rank-shape coordinate

At criticality

```text
c_*(tau)=P1/(2P0)=(1-2a)/(2a).
```

Using (1.3),

```text
boxed: c_*(tau)=cot^2[D_F(tau)/2].                        (2.1)
```

Thus all three descriptions are equivalent:

```text
maximize P0
<=> minimize c_*
<=> maximize Fisher distance from the rank-one cusp.       (2.2)
```

The #781 hexagonal conjecture can therefore be stated geometrically:

> among critical torus moduli, the hexagonal torus is the rank law farthest from the pure rank-one cusp in Fisher geometry.

No thermal metric enters this statement.

## 3. Source-zero collision is the Fisher midpoint

The bounded rank-source partition function at criticality is

```text
Z_X(s)=[c_*+cosh s]/[c_*+1].
```

The collision condition is `c_*=1`, equivalently

```text
P0=P2=1/4,
P1=1/2.
```

By (2.1),

```text
boxed: c_*=1  <=>  D_F=pi/2.                              (3.1)
```

So the modular source-zero collision locus is exactly the intersection of the critical Fisher meridian with the geodesic sphere of radius `pi/2` around the rank-one cusp.

The three zero regimes become

```text
D_F < pi/2  <=> c>1 : zeros at i*pi +- arcosh(c),
D_F = pi/2  <=> c=1 : double zero at i*pi,
D_F > pi/2  <=> c<1 : nearest zero pair purely imaginary. (3.2)
```

This gives the collision curve an intrinsic information-geometric meaning; it is still not a thermodynamic phase boundary or a Jordan diagnostic.

## 4. Numerical controls

Using the independent Arguin/Pinson controls in `critical-rank-modular-morse-20260914.md`:

### Hexagonal modulus

```text
P0 = 0.316053412701401686065271069685...
D_F = 2 asin sqrt(2P0)
    ~= 1.8381848356043984
    ~= 0.585112405806006 pi.
```

### Square modulus

```text
P0 = 0.309526275429831322769693633353...
D_F ~= 1.8112106992599288
    ~= 0.576526271536292 pi.
```

Both lie on the imaginary-zero side of the `D_F=pi/2` collision circle, with the hexagonal point farther from the rank-one cusp.

## 5. Differential consequences

Since `c=c(D_F)` is strictly decreasing,

```text
nabla c = -csc^2(D_F/2) cot(D_F/2) * nabla D_F.           (5.1)
```

At any modular stationary point of the critical rank law, `nabla c=0` iff `nabla D_F=0`. At such a point the Hessians differ only by the positive/negative scalar derivative `dc/dD_F`:

```text
Hess c = (dc/dD_F) Hess D_F.                              (5.2)
```

Therefore the Morse classification in the modular note can equivalently be read as a classification of extrema of Fisher distance.

## 6. Claim boundary

- Exact: (1.1)--(3.2), (5.1)--(5.2).
- Numerical: the displayed square/hex distances.
- Conjectural: the global hexagonal maximum of `D_F`, equivalent to the #781 global minimum of `c_*`.
