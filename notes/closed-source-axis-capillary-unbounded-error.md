# Unbounded capillary control and the first possible crossover

## Result

Let `r=m^-1`, `c=L/m`, and introduce

```text
beta = L/m^2 = c/m,
alpha = c beta = L^2/m^3.                                  (1)
```

The finite-c law from the complete two-gas transfer extends to an
unbounded capillary window under the sufficient gate

```text
m->infinity,       c may diverge,       alpha->0.            (2)
```

Uniformly in this regime,

```text
Ustar/A_N
 =-[L^2 m^(-(2L+1))/Delta]
   [I0(2c)^2-I1(2c)^2]
   [1+O(alpha+1/m+L/m^3+exp(-kappa Lm))].                   (3)
```

In particular original U remains negative.  When c also tends to infinity,

```text
Ustar/A_N
 =-m^(-(2L-1)) exp(4c)/(8 pi Delta)
   [1+O(c^-1+alpha)].                                       (4)
```

The parameter alpha is the first **rigorously possible** failure scale,
but it comes from only one mechanism: horizontal reversals which destroy
the one-height-per-column representation near the hard endpoint.  Two
other apparent threats do not saturate it.

- Arbitrarily long directed vertical runs retain the determinant and have
  relative error `O(1/m+L/m^3)`.
- Narrow-width collision is already exactly the Karlin--McGregor factor
  `J_1=I0^2-I1^2`; vertical wrapping is `exp[-Omega(Lm)]`.

At `alpha=O(1)` the bulk density of horizontal-reversal defects still
vanishes (`beta=alpha/c`), but one defect is comparable to the endpoint
survival probability `J_1/I0^2~1/(2c)`.  The new object is consequently a
single-hairpin contact kernel, not a new bulk capillary gas.  This note
defines that kernel precisely; its coefficient is not inferred from the
positive directed subfamily.

There is also a separate fixed-m result.  The exact directed positive
subfamily has exponential base

```text
q_dir(m)=[(m+1)/(m(m-1))]^2,                                 (5)
```

and proliferates when `1<m<1+sqrt(2)`.  This rules out a positive bare
axis interface tension based only on straight stripes in that range.  It
does not prove a rank-one phase or determine original U.  For
`m>1+sqrt(2)`, decay of this subfamily gives no upper bound on the complete
interface sector.

No Monte Carlo, finite-size fit or new observable enters the argument.

The logical status is deliberately split.  Equations (10)-(16), (20)-(21)
and (27)-(28) are exact-kernel/local-CLT statements.  Equation (18) is a
rigorous relaxed-word upper bound, and together with the carrier
suppression from the finite-c proof gives the sufficient full-lattice gate
`alpha->0` in (3).  Equations (23)-(25) identify the first competing
packet but do not evaluate it; whether `alpha` is sharp is therefore a
mechanism question, not part of the theorem.  If its stripped endpoint
coefficient vanishes by a Dirichlet cancellation, the gate weakens from
`alpha->0` to the bulk condition `beta->0`; no such cancellation is
assumed here.

## 1. Why an endpoint error is amplified

Write

```text
J_1(c)=I0(2c)^2-I1(2c)^2,
rho(c)=J_1(c)/I0(2c)^2.                                    (6)
```

Standard Bessel ratio bounds, or one-dimensional Laplace estimates, give
uniform positive constants `a0,A0` such that

```text
a0/(1+c) <= rho(c) <= A0/(1+c),
rho(c)=1/(2c)+O(c^-2)        as c->infinity.                 (7)
```

Thus an exceptional family of relative bulk weight epsilon has only the
crude endpoint bound `O[(1+c)epsilon]`.  Long vertical runs and horizontal
reversals each have bulk activity at most

```text
epsilon_bulk=O(beta)                                        (8)
```

by a union bound over L columns.  Dividing (8) by (7) gives

```text
O[(1+c)beta]=O(beta+alpha).                                 (9)
```

Equation (9) is the source of alpha.  It is not obtained by changing the
raw winding exponent; it is the amplification of a rare bulk defect by a
hard endpoint whose survival probability is only `O(1/c)`.

## 2. Long vertical runs preserve the determinant and are smaller

For one directed boundary the exact arbitrary-run symbol is

```text
t_r(z)=sum_(d in Z)r^|d|z^d
      =(1-r^2)/[(1-rz)(1-r/z)].                              (10)
```

Let

```text
P_d(L,r)=[z^d]t_r(z)^L.                                     (11)
```

Planarity and Lindstrom--Gessel--Viennot give the exact two-boundary
endpoint partition

```text
J_run(L,r)=P_0(L,r)^2-P_1(L,r)^2.                            (12)
```

Long runs therefore cannot add an unsuppressed `I0^2` term: they preserve
the Dirichlet zero at coincident boundaries.

For unbounded c the normalized increment law has

```text
t_r(1)=(1+r)/(1-r),
V=L Var(d)=2Lr/(1-r)^2=2c/(1-r)^2.                           (13)
```

The uniform local central limit theorem gives

```text
P_0=t_r(1)^L/[sqrt(2 pi V)] [1+O(V^-1)],
P_1/P_0=exp[-1/(2V)][1+O(V^-2)],

J_run=t_r(1)^(2L)/(2 pi V^2)[1+O(V^-1)].                    (14)
```

Comparing the two saddle expansions before truncating their common
Edgeworth terms yields

```text
J_run/J_1
 =(1-r)^4 exp[4Lr^3/3+O(Lr^5)]
   [1+O(r/c+r^2)].                                         (15)
```

and hence

```text
J_run/J_1=1+O(1/m+L/m^3).                                  (16)
```

The cancellation of the naked `O(c^-1)` term is important: it is the
same lattice local-CLT correction in numerator and denominator, not an
arbitrary-run effect.  For bounded c, the same statement follows directly
by perturbing the Fourier symbol in (10); the error is `O(Lr^2)=O(1/m)`.
Thus the compact-c and saddle regimes overlap and give (16) uniformly.

Even when `alpha=O(1)`, one has `L/m^3=alpha/L->0`.  Arbitrary directed
vertical runs are therefore not the alpha crossover.

## 3. Horizontal reversals supply the only alpha-scale bound

For one oriented contour let W=w be the number of west steps and let
U=D=n.  Its excess length is `2(w+n)`.  If simplicity is forgotten, its
based weight is bounded by

```text
C_L(w,n)r^(2w+2n),
C_L(w,n)=(L+2w+2n)!/[(L+w)!w!n!^2].                         (17)
```

For w at least one, summing (17) relative to the directed bridge gives,
uniformly while beta is small,

```text
Z_reversal,bulk/Z_directed,bulk
 <= C beta exp(C beta).                                      (18)
```

The uniform comparison behind (18) is explicit.  At fixed vertical
variation n,

```text
C_L(w,n)/C_L(0,n)
 =1/w! product_(j=1)^(2w)(L+2n+j)
       /product_(j=1)^w(L+j).                               (18a)
```

After multiplication by `r^(2w)`, this is `beta^w/w!` times
`exp[O(w(n/L+w/L))]`.  Under the `w=0` bridge weight, `n=O(c)` with an
exponential tail, while `c/L=1/m`; summing first over n and then over w
gives `exp(C beta)-1`.  The complementary region `n+w>L/4` is controlled
by the same multinomial factorial tail and is `exp[-Omega(L log m)]`.
Thus (18) does not assume `L^2` independent planar positions: a west step
has only one uncompensated column coordinate, so its activity is
`Lr^2=beta`, not `L^2r^2`.

In particular the first competing packet is not schematic.  Its exact
fixed-n multiplier over the directed bridge is

```text
r^2 C_L(1,n)/C_L(0,n)
 =r^2 (L+2n+1)(L+2n+2)/(L+1)
 =beta[1+O(n/L+1/L)].                                      (18b)
```

Thus `w=1` is the unique bulk packet capable of reaching the amplified
`alpha` scale; every `w>=2` packet is `o(c^-1)` when alpha is bounded.

The same estimate holds for either boundary and after the tight two-gas
decorations are restored.  Combining (18) with (7) proves

```text
|Z_endpoint,full-J_run|/J_1
 <= C(1+c)beta exp(C beta)+exp[-kappa L^2/c]
 =O(alpha+beta)+exp[-kappa Lm].                              (19)
```

This is a sufficient bound rather than a claim that every reversal
saturates it.  A physical reversal remains noncrossing, so a sharper
Dirichlet-preserving polymerization may cancel the factor c.  What is
proved here is the precise first scale at which the current encoding can
no longer discard the reversal sector while retaining an endpoint-marked
statement.

Equations (16) and (19), together with the two-gas cancellation and the
suppression of extra essential components, prove (3).

## 4. Narrow widths introduce no additional parameter

In the unit-jump limit the one-particle kernel is

```text
p_c(a,b)=I_(b-a)(2c).
```

For initial gap d, strict nonintersection is exactly

```text
J_d(c)=det[[I0(2c),I_d(2c)],
           [I_d(2c),I0(2c)]].                               (20)
```

Thus all repeated near-collisions at width one have already been summed
inside J_1.  They are not a perturbative `O(1/L)` event after the endpoint
occupation leverage is inserted.  The only finite-L alias is winding of
the height bridge around the transverse circle.  Its Gaussian tail is

```text
O[exp(-kappa L^2/c)]=O[exp(-kappa Lm)].                      (21)
```

This proves the last error term in (3).

## 5. The alpha=O(1) crossover is a single-hairpin contact kernel

Suppose `c->infinity`, `beta->0`, but `alpha=c beta` has a finite nonzero
limit.  Two reversal defects have bulk weight `O(beta^2)` and endpoint
relative weight

```text
c beta^2=alpha^2/c->0.                                      (22)
```

Only the one-defect sector can survive.  Root a reversal at its leftmost
west edge and cut the contour immediately before and after its irreducible
hairpin.  Let

```text
H_c^phys(d_in,d_out)
```

be the resulting positive two-boundary transfer: unit-jump noncrossing
propagation before and after the rooted hairpin, summed over its finite
vertical span and over which boundary carries it.  The physical
occupied-corner resolution is imposed inside H, so crossing and a shared
cut edge have zero weight.  Here `H^phys` includes the physical reversal
activity.  Strip its single bulk factor by setting

```text
Hhat_c=beta^(-1) H_c^phys.
```

Normalize the stripped kernel by the bulk `I0(2c)^2` factor and define
the endpoint coefficient

```text
kappa_hp(c)
 =I0(2c)^(-2)<d=1| Hhat_c |d=1>
 =beta^(-1) I0(2c)^(-2)
   <d=1| H_c^phys |d=1>.                                   (23)
```

Then the exact first possible crossover form is

```text
Z_endpoint/I0(2c)^2
 =rho(c)+beta kappa_hp(c)+o(c^-1),                           (24)
```

or, if `kappa_hp(c)` has a finite limit,

```text
c Z_endpoint/I0(2c)^2
 ->1/2+alpha kappa_hp(infinity).                             (25)
```

Equations (23)-(25) are a concrete new kernel, not a fitted scaling
function.  Its unresolved datum is the signed/positive coefficient of the
one irreducible physical hairpin.  A proof that the hard noncrossing
Dirichlet zero forces `kappa_hp(infinity)=0` would show alpha is only a
crude gate and extend (3) to beta->0.  A nonzero value would give the first
genuine unbounded-capillary crossover.  No bulk resummation beyond one
defect is needed at fixed alpha.

## 6. Fixed-m directed roughening

The same exact one-boundary kernel also answers a different limit order.
Fix `m>1` and let L grow.  The Fourier maximum of (10) is

```text
t_r(1)=(m+1)/(m-1).
```

After restoring the straight horizontal factor `m^-L`, one boundary has
growth base `(m+1)/[m(m-1)]`; two boundaries have (5).  Therefore

```text
q_dir(m)>1  iff  1<m<1+sqrt(2),
q_dir(m)=1  iff  m=1+sqrt(2).                                (26)
```

The local central limit theorem (14), now at fixed m, gives the explicit
endpoint determinant

```text
K_0(L,m)^2-K_1(L,m)^2
 =[(m+1)/(m-1)]^(2L)
   (m-1)^4/[8 pi m^2 L^2]
   [1+O(L^-1)].                                             (27)
```

Including the closed-source straight cost `m^(-(2L-1))` turns this into

```text
m^(-(2L-1))[K_0^2-K_1^2]
 =(m-1)^4/[8 pi m L^2] q_dir(m)^L[1+O(L^-1)].               (28)
```

For `m<1+sqrt(2)`, (28) is a rigorous positive subfamily which grows
exponentially with L.  It excludes:

- a Peierls argument that uses the straight cost `m^(-2L)` with only a
  subexponential degeneracy;
- a claim that the directed rank-one interface has positive bare tension.

It does **not** prove a rank-one thermodynamic phase, a nonzero rank-one
probability, or a sign of original U: the rank-zero/rank-two bulk
partitions have area-order free energies and normalization still matters.

For `m>1+sqrt(2)`, (28) decays.  Since it is a lower bound supplied by one
positive subfamily, that decay cannot exclude proliferation from horizontal
reversals, multiple carriers or another interface class.  The threshold is
therefore a sharp directed-subfamily statement, not an upper bound for the
full fixed-m model.

## Scientific card

- **Full-law gate:** `alpha=L^2/m^3->0` is sufficient for the Bessel
  endpoint law and negative original U even when `L/m->infinity`.
- **Mechanism isolation:** long vertical runs and repeated narrow
  collisions preserve the determinant; only horizontal reversals can
  currently saturate alpha.
- **New crossover object:** the one-hairpin endpoint kernel (23)-(25).
  It decides whether alpha is sharp or merely conservative.
- **Independent fixed-m result:** the directed positive subfamily changes
  exponential sign at `m=1+sqrt(2)` with prefactor (27)-(28).
- **Boundary:** no claim is made for `alpha=O(1)` until the one-hairpin
  coefficient is evaluated, or for the full fixed-m interface from a
  directed lower bound alone.
