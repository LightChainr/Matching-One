# The delocalized one-west word has no endpoint pole

## Result

The relaxed one-west packet can borrow vertical edges from the entire rough
bridge, so it is not the finite local operator studied in
`closed-source-single-hairpin-endpoint-zero.md`. Nevertheless its exact
word generating function shows that the induced gap kernel is short-range.

After stripping the declared bulk activity

```text
beta=Lr^2=L/m^2,       r=m^-1,                                (1)
```

the one-west kernel is a common analytic Toeplitz multiplier

```text
b_r(z)=(1+2/L)[1-r(z+z^-1)]^-2.                              (2)
```

It has neither

```text
1/[(p_2-p_1)(q_2-q_1)]
```

nor any weaker low-momentum pole capable of cancelling the two external
antisymmetric zeros. In the joint limit `r->0`, `c=Lr->infinity`,
`beta->0`, its endpoint coefficient obeys

```text
kappa_w(c)=c^-1[1+o(1)],
kappa_w(infinity)=0,
c kappa_w(c)->1.                                              (3)
```

Thus even the fully delocalized relaxed one-west word changes the directed
endpoint only by relative order beta. It cannot realize an alpha-scale
packet with `alpha=c beta`. Borrowing a long interval of bridge history is
not the same as producing a long-range gap operator.

No physical self-avoidance assumption is needed for this conclusion: it is
already true in the larger relaxed-word family.

## 1. Exact no-west and one-west word kernels

Fix net horizontal displacement L and mark vertical displacement with z.
With no west step, a word has L east steps, u up steps and d down steps.
Summing every ordering gives

```text
D_L(z)
 =sum_(u,d>=0) (L+u+d)!/[L!u!d!]
                r^(u+d) z^(u-d)
 =[1-r(z+z^-1)]^(-(L+1)).                                    (4)
```

Write

```text
g_r(z)=1-r(z+z^-1).
```

With exactly one west step there must be `L+1` east steps. The exact relaxed
word sum is

```text
W_L(z)
 =sum_(u,d>=0) (L+u+d+2)!/[(L+1)!u!d!]
                r^(u+d+2) z^(u-d)
 =r^2(L+2) g_r(z)^(-(L+3)).                                  (5)
```

For zero net vertical displacement, the coefficient ratio at fixed
`u=d=n` is consequently

```text
r^2 C_L(1,n)/C_L(0,n)
 =r^2 (L+2n+1)(L+2n+2)/(L+1),                                (6)
```

recovering the relaxed contour-word multiplier exactly.

Divide (5) by beta. Equations (1), (4) and (5) give

```text
What_L(z):=beta^-1 W_L(z)
 =(1+2/L)g_r(z)^(-2)D_L(z),                                  (7)
```

which proves (2). The two extra powers of `g_r^-1` encode every possible
borrowing of vertical edges. They do not create a resolvent in the gap.

## 2. The apparent long word is a short-range gap kernel

Let

```text
s=sqrt(1-4r^2),
q=(1-s)/(2r)=r+O(r^3).                                       (8)
```

The Fourier coefficients of the multiplier in (2), without the harmless
factor `1+2/L`, are exact:

```text
[z^ell]g_r(z)^(-2)
 =q^|ell|[(|ell|+1)-(|ell|-1)q^2]
   /[s^2(1-q^2)].                                             (9)
```

Hence its gap tail is exponential, with correlation length
`1/|log q|`, and

```text
sum_ell |[z^ell]g_r(z)^(-2)|
 =g_r(1)^(-2)=(1-2r)^(-2).                                  (10)
```

In particular the stripped gap kernel is uniformly summable as `r->0`.
The west step is delocalized only in its choice of word positions; after
those positions are summed, its action on the gap is more local, not less.

## 3. Antisymmetric Fourier kernel

Let

```text
A_n=[z^n]D_L(z),       B_n=[z^n]What_L(z).
```

Putting the west step on either boundary and antisymmetrizing gives the
endpoint matrix element

```text
Hhat_w(1)=2[A_0B_0-A_1B_1].                                  (11)
```

Equivalently,

```text
Hhat_w(1)
 ={1/(2pi)^2} integral integral
   D_L(e^(ip))D_L(e^(iq))
   [b_r(e^(ip))+b_r(e^(iq))]
   [1-cos(p-q)] dp dq.                                       (12)
```

The last factor is the squared external Vandermonde:

```text
1-cos(p-q)=(p-q)^2/2+O((p-q)^4).                             (13)
```

At the joint saddle `p,q=O(c^-1/2)`,

```text
b_r(e^(ip))
 =(1+2/L)(1-2r)^(-2)[1+O(rp^2)].                             (14)
```

It is finite and even. Thus (12) retains the full factor (13). To produce
an unrestricted-bulk contribution, a two-sided insertion would need a
kernel singular as

```text
[(p_2-p_1)(q_2-q_1)]^-1,                                    (15)
```

or, after the diagonal reductions leading to (12), as `(p-q)^-2`.
Equations (9) and (14) exclude both. There is no fractional or logarithmic
singularity either, because (2) is analytic in an annulus containing the
unit circle for all sufficiently small r.

## 4. Scaling coefficient

The normalized no-west word distribution associated with (4) has variance

```text
V_L=2(L+1)r/(1-2r)=2c[1+O(r+L^-1)].                          (16)
```

The cancellation-sensitive local CLT gives

```text
A_0^2-A_1^2
 =D_L(1)^2/[2pi V_L^2][1+O(V_L^-1)].                         (17)
```

Across the saddle, (14) is constant to relative `o(1)`, so (11) becomes

```text
Hhat_w(1)
 =2(1+2/L)(1-2r)^(-2)
   [A_0^2-A_1^2][1+o(1)].                                   (18)
```

When `r->0`, `c->infinity` and `beta=Lr^2->0`,

```text
D_L(1)^2=exp(4c)[1+O(beta+r)],
V_L=2c[1+O(r+L^-1)],
I0(2c)^2=exp(4c)/(4pi c)[1+O(c^-1)].                         (19)
```

Combining (17)-(19) yields

```text
I0(2c)^(-2) Hhat_w(1)
 =c^-1[1+O(c^-1+beta+r)],                                    (20)
```

which is (3). The factor one in `c kappa_w(c)->1` includes both carrier
boundaries; one carrier alone contributes one half.

## 5. Consequence for the proposed alpha crossover

Restoring the one-west activity beta gives

```text
beta kappa_w(c)=beta/c[1+o(1)].                               (21)
```

Since the directed endpoint is

```text
rho(c)=1/(2c)+O(c^-2),                                       (22)
```

the relaxed one-west response has the finite relative size

```text
[beta kappa_w(c)]/rho(c)=2beta[1+o(1)].                       (23)
```

It therefore vanishes whenever `beta=L/m^2->0`, including the putative
`alpha=O(1)`, `c->infinity` crossover where `beta=alpha/c->0`.

The strict zero criterion is now explicit:

> A one-west family can contribute at alpha scale only if, after stripping
> beta and summing its borrowed bridge segments, its antisymmetric Fourier
> kernel develops an inverse external-Vandermonde pole. The exact relaxed
> word kernel (2) is analytic and fails this criterion.

Any surviving alpha mechanism must therefore use extra nonlocal structure
not contained in one west step plus arbitrary borrowed vertical edges—for
example a jointly scaled collar state whose gap range itself grows like
`sqrt(c)`, or an operation which opens the noncrossing wall. Neither is
present in the relaxed one-west word.

## Scientific card

- **Exact relaxed kernel:** equations (4)-(7) sum all placements of one west
  step and all borrowed vertical bridge edges.
- **Gap diagnosis:** its additional kernel is the exponentially decaying
  coefficient sequence (9), not a long-range resolvent.
- **Low-momentum result:** the multiplier is analytic and leaves the squared
  Vandermonde (13) intact; no alpha-producing pole occurs.
- **Scaling:** `kappa_w(c)~1/c`, with `c kappa_w(c)->1` after summing both
  carrier boundaries.
- **Boundary:** this is the relaxed one-west envelope. A genuinely nonlocal
  collar state with range growing as `sqrt(c)` is a different operator and
  is not ruled out here.
