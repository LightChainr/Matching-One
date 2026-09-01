# The complete axis one-carrier reversal sector has a beta gate

## Result

Let `Z_1c(L,m)` be the axis endpoint weight of one essential carrier with:

- any number of horizontal reversals;
- physical self-avoidance imposed by deleting illegal words;
- the matched black/white local clouds;
- a finite set of collar states of fixed lattice radius R;
- crossings and a shared cut edge forbidden.

Let `J_run(L,m)` be the exact directed noncrossing endpoint determinant with
arbitrary vertical runs. Under the minimal locality conditions stated in
Section 4,

```text
|Z_1c-J_run|/J_run <=C_R beta,       beta=L/m^2,              (1)
```

uniformly in the capillary parameter `c=L/m`. Consequently

```text
beta->0       implies       Z_1c/J_run ->1.                   (2)
```

Together with the exact-run comparison, this gives

```text
Z_1c/J_Bessel
 =1+O(beta+m^-1+L/m^3),                                      (3)
```

where `J_Bessel=I0(2c)^2-I1(2c)^2`. Since `L/m^3=beta/m`, every error in
(3) vanishes under `beta->0`.

Thus the complete horizontal-reversal part of the axis one-carrier gate is

```text
L/m^2 ->0,                                                    (4)
```

not `L^2/m^3->0`. The stronger alpha gate arose from comparing an arbitrary
bad-word mass with the unrestricted bulk `I0^2` and forgetting that physical
deletion preserves the two external Dirichlet zeros.

This result is deliberately one-carrier. It does not control an additional
essential component or a collar radius that grows with c.

## 1. The all-word analytic majorant

The exact relaxed word with w west steps is

```text
W_(L,w)(z)
 =r^(2w) binom(L+2w,w)
  [1-r(z+z^-1)]^(-(L+2w+1)),       r=m^-1.                   (5)
```

Summing every `w>=0` gives the no-west kernel times the analytic multiplier

```text
M_L(x)=[1-4x]^-1/2
       {2/[1+sqrt(1-4x)]}^L,
x=r^2/[1-r(z+z^-1)]^2.                                      (6)
```

For `beta=Lr^2->0`, there is an annulus about the unit circle, independent
of L, on which

```text
M_L(x)=1+O(beta),                                             (7)
```

and all weighted coefficient moments of the **relative multiplier** are
bounded by `C beta`. In particular, if

```text
e_L(n)=[z^n]{M_L(x(z))-1},
```

then for a fixed small `eta>0`,

```text
sum_n exp(eta|n|)|e_L(n)| <=C_eta beta.                       (8)
```

The full relaxed reversal transfer is the directed kernel convolved with
`e_L`. Equation (8) is stronger than a total-mass estimate for the inserted
packet: it rules out a relative gap tail whose range grows on the capillary
scale.

## 2. Why deleting illegal words cannot create a pole

Erase the reversal record from a physical word, leaving its directed
backbone. The contour-word rooting used in the all-w count is injective after
the backbone, reversal number and cut data are retained. Physical
self-avoidance then only deletes records from this positive relaxed packet
majorant; it does not assign a larger weight to a survivor. Thus the packet
kernel relative to each directed backbone is coefficientwise bounded by
(8). Antisymmetric signs appear only after the two labelled carriers are
combined, so the absolute value remains bounded by the sum of the positive
labelled majorants.

There are two geometrically distinct pieces.

### 2.1 Away from the other boundary

If the carrier packet and its radius-R collar do not meet the other boundary,
vertical translation of the whole word is exact. Self-avoidance may delete a
complicated global subset, but the surviving one-particle kernel remains a
Toeplitz convolution. Coefficientwise domination by (8) makes its Fourier
symbol analytic and

```text
symbol_physical/symbol_directed=1+O(beta).                    (9)
```

The same multiplier enters both rows of the two-path determinant. Hence its
endpoint effect is coherent and relative `O(beta)`; it is not an additive
`O(beta I0^2)` error.

### 2.2 Near the hard wall

Dependence on the spectator boundary is possible only when a packet/collar
reaches the gap wall. Fixed R and the geometric vertical-span tail imply an
exponentially localized two-gap kernel `B_wall` with

```text
sum_(d,e>=1) exp[eta(d+e)] |B_wall(d,e)| <=C_R beta.          (10)
```

Deleting additional illegal words only decreases the positive majorant in
(10). Its double sine transform therefore satisfies

```text
Bhat_wall(p,q)=O(beta p q)       as p,q->0.                   (11)
```

The factors p and q are the incoming and outgoing Vandermonde zeros. An
inverse factor `1/(pq)` would require divergent first gap moments, which
contradicts (10). Thus deletion can make the wall kernel non-Toeplitz, but it
cannot manufacture an inverse-Vandermonde pole.

This is the exact answer to the deletion concern: a subset of a positive
exponentially summable family remains exponentially summable. Removing words
cannot turn an analytic generating function into a threshold resolvent.

## 3. Endpoint comparison

The directed endpoint obeys the uniform bound

```text
a/(1+c) <=J_run/bulk_weight<=A/(1+c).                         (12)
```

The common Toeplitz part of Section 2.1 changes both numerator and bulk
weight by `1+O(beta)`. For the wall part, equations (10)-(11) and the
Dirichlet saddle give

```text
|delta Z_wall|<=C_R beta bulk_weight/(1+c).                   (13)
```

Dividing (13) by (12), and adding the coherent bulk correction, proves (1).
No factor c appears because both terms retain the same endpoint zero as
`J_run`.

The proof also covers any finite number of internal collar labels. Write
their transfer as a fixed-dimensional block matrix. Each block satisfies
(8) or (10), and a finite matrix sum changes only `C_R`, not the low-momentum
power.

## 4. Minimal conditions and failure modes

The conclusion requires exactly the following conditions.

1. **Backbone domination.** Erasing reversals maps each physical word
   injectively to a directed backbone plus a relaxed reversal record;
   physical filtering only deletes such records, or multiplies them by
   uniformly bounded nonnegative local cloud weights.
2. **Wall preservation.** Crossing and a shared cut edge remain illegal, so
   the ordered transfer never opens the `d=0` Dirichlet wall.
3. **Fixed-radius collar.** The number of collar labels and their radius R do
   not grow with L,m or c; extra vertical span retains its geometric r-tail.
4. **One carrier.** No additional essential component supplies a second
   delocalized coordinate.

These are minimal in the following operational sense. If condition 1 fails,
signed reweighting can cancel the analytic majorant. If condition 2 fails,
the constant bulk mode becomes accessible. If condition 3 fails with
`R~sqrt(c)`, the first moments in (10) need not be uniform and a threshold
kernel may emerge. If condition 4 fails, the one-gap reflection reduction no
longer describes the state space.

Global self-avoidance by itself is not an extra condition: far from the wall
it preserves translation covariance, and near the wall it only deletes terms
already bounded by (10).

## 5. Scientific boundary

Equation (1) closes the axis one-carrier horizontal-reversal problem under a
fixed-radius physical collar. It does not close:

- a collar or marked state whose range grows like `sqrt(c)`;
- two or more essential carriers;
- a nonlocal black/white-cloud interaction not dominated by (5);
- the separate tilted-sector and sector-odds assumptions entering original
  U.

Subject to those independent boundaries, horizontal reversals no longer
justify the alpha gate. Their complete contribution is endpoint-relative
`O(beta)=O(L/m^2)`.

## Scientific card

- **All reversals:** the Catalan multiplier (6) supplies a positive analytic
  majorant with total excess `O(beta)`.
- **Self-avoidance:** coefficientwise deletion preserves exponential gap
  summability and cannot create an inverse-Vandermonde pole.
- **Collar closure:** every fixed-radius wall correction retains both sine
  zeros and contributes `O(beta/(1+c))` in bulk normalization.
- **New gate:** the complete axis one-carrier reversal sector needs only
  `L/m^2->0` for relative endpoint error `o(1)`.
- **Remaining boundary:** growing-range collars and extra essential
  components are not covered.
