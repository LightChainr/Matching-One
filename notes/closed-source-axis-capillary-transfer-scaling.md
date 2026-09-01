# Matched two-cloud capillary transfer gives a positive Bessel determinant

**Promotion.**  This note derives the determinant inside the matched
directed model and records the correction chain.  The companion
[complete transfer](closed-source-axis-capillary-complete-transfer.md) and
[signed interface theorem](closed-source-axis-signed-interface-transfer.md)
supply the overhang, endpoint and extra-component estimates that promote
the same formula to the leading square-lattice axis law for bounded
`L/m`.

## New result and correction history

In the capillary window

```text
L,m -> infinity,       c=L/m in [0,infinity),
a=m^-2=c^2/L^2,        h_*=1+a,                                (1)
```

the directed rank-one calculation must retain both local clouds:

1. occupied black singletons outside the essential stripe;
2. vacant white single holes inside it.

After matching these two clouds at `h=h_*`, and then imposing noncrossing
of the two directed stripe boundaries, the scaling coefficient is

```text
Phi_2cloud(c)
 =J_1(c)/Delta
 =[I0(2c)^2-I1(2c)^2]/Delta.                                  (2)
```

Thus the corresponding original response is

```text
Ustar/A_N
 =-L^2 m^(-(2L+1)) Phi_2cloud(c)[1+o(1)].                     (3)
```

Because `I0(2c)>I1(2c)>=0`, this matched two-cloud directed candidate is
strictly positive for every finite `c>=0`; hence (3) stays negative and has
no capillary zero in this model.

This supersedes two earlier incomplete candidates in this note's history:

- multiplying the rigid result only by `I0(2c)^2` missed all derivative
  leverage of the width mode;
- retaining only the exterior black-singleton cloud produced the restricted
  bracket

```text
1-c^2/6+(2c/3) I1(2c)/I0(2c),                                (4)
```

  and an apparent zero. Equation (4) is withdrawn as a physical two-cloud
  result. Its `-c^2/6` term is an artefact of omitting the interior white-hole
  pressure, while its separate bridge-area term is cancelled by the matched
  two-cloud width tilt.

The value of the correction is mechanistic: the final positivity is not the
naive positivity of a partition multiplier. It follows from an exact local
cloud match followed by a noncrossing two-boundary determinant.

No simulation, finite-c fit, new source or lattice enumeration is used.

## 1. One directed boundary

Represent an axis-essential boundary by its integer height after each of L
forward columns. With `r=m^-1`, its exact displacement kernel is

```text
T_r(x,y)=r^|x-y|.                                               (5)
```

The Fourier eigenvalue is

```text
t_r(theta)=sum_(d in Z) r^|d| exp(i d theta)
          =(1-r^2)/(1-2r cos(theta)+r^2).                      (6)
```

After fixing the global height, the closed-boundary partition is

```text
Z_1,L(r)
 =[z^0](sum_(d in Z) r^|d| z^d)^L
 =(1/2pi) integral_0^(2pi) t_r(theta)^L dtheta.                 (7)
```

For `r=c/L`, uniformly at fixed c,

```text
L log t_(c/L)(theta)=2c cos(theta)+O(c^2/L),
Z_1,L(c/L) -> I0(2c).                                          (8)
```

More generally, a directed bridge with net height displacement d has
scaling kernel `I_d(2c)`. This off-diagonal entry is essential at narrow
stripe width; retaining only the diagonal factor `I0(2c)` treats the two
boundaries as independent even when they are one lattice spacing apart.

## 2. The omitted local cloud

For a straight stripe of width w, define the two available local-cloud
areas

```text
H_w=L max(w-2,0),
M_w=L max(L-w-2,0).                                            (9)
```

`H_w` counts positions of the interior white-hole cloud and `M_w` those of
the exterior black-singleton cloud. Their combined stripe weight is

```text
W_w(h)
 =h^(Lw) (1+a^2/h)^(H_w) (1+a h)^(M_w).                       (10)
```

Introduce

```text
C(h)=h+a^2,             B(h)=1+a h.                            (11)
```

At the matched root,

```text
C(h_*)=h_*+a^2=1+a+a^2=B(h_*).                                (12)
```

For every bulk width `2<=w<=L-2`, equations (9)-(12) give

```text
W_w(h_*)
 =h_*^(2L) B(h_*)^[L(L-4)],                                   (13)
```

which is independent of w. Equivalently, the width fugacity is exactly

```text
q_width=C(h_*)/B(h_*)=1.                                      (14)
```

This equality is the missing cancellation. The exterior-only calculation
used

```text
h^(Lw) B(h)^[L max(L-w-2,0)]
```

and obtained the normalized logarithmic response

```text
c^2/2-c^4/12.                                                  (15)
```

The `c^2/2` came from the width endpoint, while `-c^4/12` was the unmatched
bulk pressure. Restoring the interior factor in (10) supplies the opposite
bulk pressure. Directly differentiating the three width ranges
`w=1`, `2<=w<=L-2`, `w=L-1` gives instead

```text
d_h log[Z_1/(Z_0+Z_2)] -> c^2/2.                              (16)
```

In the full root/normalization combination, the width q-tilt and the
bridge-area tilt cancel. Thus neither the exterior-only `-c^4/12` nor the
previous additive bridge variance survives as an independent bulk response.
With the fixed-L prefactors, (16) restores the rigid endpoint coefficient;
there is no physical `1-c^2/6` factor in the matched two-cloud model.

This also explains why an endpoint of only `O(1/L)` width weight cannot be
discarded: its occupation leverage is `O(L^2)`. The lesson from the previous
restricted model remains valid, but its signed coefficient does not.

## 3. Noncrossing resummation at the endpoint

At a bulk width, two boundary transfers factor to `I0(2c)^2`. At the
minimal allowed separation d, however, the two directed boundaries must not
cross. The two-path transfer is the determinant

```text
J_d(c)
 =det [[I0(2c), I_d(2c)],
       [I_d(2c), I0(2c)]]
 =I0(2c)^2-I_d(2c)^2.                                        (17)
```

The physical width endpoint has `d=1`. Therefore noncrossing multiplies the
straight endpoint coefficient by

```text
rho(c)=J_1(c)/I0(2c)^2
      =1-[I1(2c)/I0(2c)]^2.                                  (18)
```

The bulk two-cloud response and bridge-area terms have already cancelled in
Section 2. The surviving rigid endpoint coefficient is dressed by the bulk
partition `I0(2c)^2` and attenuated by (18). Their product is exactly

```text
I0(2c)^2 rho(c)=J_1(c),                                       (19)
```

which proves (2) within the matched, noncrossing directed model.

Equation (17) is the missing narrow-width interaction in the earlier
restricted calculation. It is not an arbitrary signed correction: it is the
minimal two-boundary exclusion determinant.

## 4. Root and denominator propagation

Let `P1(h)` be the separately normalized axis rank-one probability. The
fixed-L reciprocal cancellation and the two-cloud root give, before the
capillary resummation,

```text
P1,axis,h
 =+(L^4/2)m^(-(2L+1))[1+o(1)].                               (20)
```

The matched capillary replacement is therefore

```text
P1,axis,h
 =+(L^4/2)m^(-(2L+1)) J_1(c)[1+o(1)].                        (21)
```

The tilted rank-one sector remains exponentially later under the declared
sector-odds assumption, and the within-geometry denominator is

```text
D_h=L^2/2[1+o(1)]>0.                                         (22)
```

Finally `E=1-P1` supplies the minus sign, and the positive geometry factor
`Delta` gives (2)-(3). Root motion and denominator normalization do not
change the winding exponent or the determinant sign; their bulk capillary
pieces are already included in the two-cloud cancellation of Section 2.

## 5. Sign and asymptotics

For `c>0`, the integral representations give

```text
I0(2c)-I1(2c)
 =(1/pi) integral_0^pi exp(2c cos theta)(1-cos theta)dtheta >0,
I0(2c)+I1(2c)>0.                                              (23)
```

Hence

```text
J_1(c)=[I0(2c)-I1(2c)][I0(2c)+I1(2c)]>0,                     (24)
```

with `J_1(0)=1`. The small-c expansion is

```text
Delta Phi_2cloud(c)
 =I0(2c)^2-I1(2c)^2
 =1+c^2+c^4/2+5c^6/36+O(c^8).                                (25)
```

At large c,

```text
Delta Phi_2cloud(c)
 =exp(4c)/(8pi c^2)[1+O(c^-1)].                               (26)
```

Thus capillary entropy exponentially dresses the amplitude but does not
create a zero. Its log growth `4L/m+O(log(L/m))` remains negligible beside
the topological tension `2L log m`, so the rigid logarithmic winding
exponent is unchanged.

## 6. Exact content and promotion boundary

**Complete in the matched noncrossing directed model:**

- every directed integer-height displacement and exact bridge closure;
- both exterior black-singleton and interior white-hole local clouds;
- the exact match `C(h_*)=B(h_*)` and the resulting bulk width/area
  cancellation;
- the straight-width endpoint response (16);
- the d=1 noncrossing two-boundary determinant (17);
- propagation through the original E sign and positive denominator.

**Required if this note is read alone:**

1. uniform suppression or exact resummation of horizontal reversals and
   overhangs;
2. control of nonlocal interactions between capillary boundaries and local
   clouds beyond the independent factors in (10);
3. control of additional essential components and contractible contour
   decorations;
4. suppression of the tilted sector at the same joint limit;
5. uniform alignment of the restricted sector odds used in (22).

The two companion proofs linked at the top control these items on compact
`L/m` intervals.  Consequently positivity of (2) rules out a capillary zero
in the leading axis law there.  The remaining asymptotic boundary is
unbounded `L/m`, not an omitted exterior-only endpoint effect.

## Scientific card

- **Corrected object:** the physical local weight is (10), not an
  exterior-singleton gas alone.
- **Exact cancellation:** `C(h_*)=B(h_*)` removes the unmatched bulk
  `-c^4/12` term and the separate bridge-area response.
- **Endpoint interaction:** noncrossing of the two directed boundaries gives
  `J_1=I0^2-I1^2`, rather than the independent product `I0^2`.
- **Signed prediction:** `Phi_2cloud=J_1/Delta` is strictly positive, so the
  matched directed model predicts negative original U with no finite-c zero.
- **Withdrawn result:** the earlier bracket (4) was a controlled
  exterior-only model, not a complete local two-cloud calculation.
- **Next boundary:** unbounded `L/m`; the compact-capillary leading axis
  law is closed by the companion estimates.
