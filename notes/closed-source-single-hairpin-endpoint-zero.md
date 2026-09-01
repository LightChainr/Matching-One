# A physical single hairpin cannot lift the Dirichlet endpoint zero

## Result

Use the one-hairpin normalization introduced in
`closed-source-axis-capillary-unbounded-error.md` at commit `1163bebe`:

```text
beta=L/m^2,
Hhat_c=beta^(-1) H_c^phys,
kappa_hp(c)=I0(2c)^(-2)<d=1|Hhat_c|d=1>.                     (1)
```

Once the local operator is fixed to be the physical rooted hairpin—one
horizontal reversal, both possible carrier boundaries summed, crossings and
a shared cut edge assigned zero—the exact reflection representation gives

```text
kappa_hp(c)=O(c^-1),
kappa_hp(infinity)=0.                                        (2)
```

Therefore the proposed first-crossover expansion sharpens to

```text
Z_endpoint/I0(2c)^2
 =rho(c)+O(beta/c)+o(c^-1),
rho(c)=1/(2c)+O(c^-2).                                       (3)
```

For this packet, the correction relative to the directed endpoint is
`O(beta)`, not `O(alpha=c beta)`. A nonzero finite limit in (1) would require
an operator which temporarily removes the absorbing wall—equivalently a
crossing/shared-edge state or an inverse-Vandermonde singularity. Those are
excluded by the definition of `H_c^phys`.

This closes the single-hairpin coefficient. It does not by itself prove
that every multi-packet or nonlocal full-lattice remainder is controlled
under only `beta->0`.

## 1. Exact noncrossing state space

The exact one-boundary column transfer is

```text
T_r(a,b)=r^|a-b|,       r=m^-1.                               (4)
```

On two labelled heights, let `R` exchange the two coordinates and put

```text
P_-=(1-R)/2,
Q_r=P_-(T_r tensor T_r)P_-.                                  (5)
```

Restricting (5) to ordered pairs `x_1<x_2` is exactly the physical
noncrossing transfer. In the gap coordinate `d=x_2-x_1`, it is the killed
half-line transfer with Dirichlet wall at `d=0`. The endpoint vector is the
antisymmetric state

```text
|partial>=(|0,1>-|1,0>)/sqrt(2).                              (6)
```

The usual reflection determinant is simply

```text
<partial|Q_r^L|partial>=K_0(L,r)^2-K_1(L,r)^2.                (7)
```

Thus the factor `rho(c)~1/(2c)` is the spectral cost of the two exterior
antisymmetric zeros at the coalesced saddle. It is not a small combinatorial
population which a local physical packet may discard.

## 2. The local hairpin operator and its normalization

Root a contour word at its unique west edge. Cut immediately before the
first edge of the irreducible packet and immediately after its last edge.
Define `B_r^phys` by the following exhaustive local rules:

1. the carrier has exactly one west edge and no other horizontal reversal;
2. all finite vertical spans are summed with their original powers of r;
3. the spectator boundary is kept and the two carrier choices are summed;
4. the occupied-corner black/white collar factors inside the packet are
   included;
5. a crossing or a shared cut edge has matrix element zero.

One reversal has excess horizontal length two, so write

```text
Btilde_r=r^-2 B_r^phys.                                       (8)
```

Additional span has a geometric r-tail. Hence `Btilde_r` has finite range in
the `r->0` limit and a uniformly summable gap kernel for all sufficiently
small r. Rules 3 and 5 give its exact reflection property

```text
Btilde_r=P_- Btilde_r P_-.                                    (9)
```

In particular, (9) is not an assumption that a generic overhang is
Toeplitz. The operator may depend arbitrarily on the finite gap and collar
state. It only says that the physical packet never opens the killed wall.

Summing the root column and stripping its bulk activity gives the precise
finite-L version of (1):

```text
Hhat_(L,r)
 =(1/L) sum_(j=0)^(L-1)
   Q_r^(L-1-j) Btilde_r Q_r^j,

H_(L,r)^phys=beta Hhat_(L,r),       beta=Lr^2.                (10)
```

The factor `1/L` in (10) matters. `beta` already contains the L possible
root columns; after division by beta the operator is a placement average,
not an unnormalized time integral.

## 3. Reflection leaves two external zeros

Let `(p_1,p_2)` and `(q_1,q_2)` be the incoming and outgoing Fourier
variables. The endpoint wavefunction from (6) is

```text
psi_partial(p_1,p_2)
 =exp(i p_2)-exp(i p_1)
 =i(p_2-p_1)+O(|p|^2).                                       (11)
```

There is one such factor on each side of (10). The dominant capillary saddle
has all momenta of order `c^-1/2`; their product supplies an unavoidable
factor `O(c^-1)` relative to the unrestricted bulk saddle `I0(2c)^2`.

The placement average cannot cancel these zeros. In a spectral pair with
transfer eigenvalues lambda and mu, its exact multiplier is the divided
difference

```text
R_L(lambda,mu)
 =(1/L) sum_(j=0)^(L-1) lambda^(L-1-j) mu^j
 =(lambda^L-mu^L)/[L(lambda-mu)],                             (12)
```

continued at `lambda=mu` by `lambda^(L-1)`. Hence (12) has no resolvent pole
at the coalesced saddle. The Fourier kernel of `Btilde_r` is bounded there
by the geometric-span estimate following (8). Combining (11)-(12) with the
same saddle bound used for the exact directed determinant gives

```text
|<partial|Hhat_(L,r)|partial>|
 <=C I0(2c)^2/(1+c)                                          (13)
```

uniformly as `r->0`, including unbounded c in the single-defect window.
Equations (1) and (13) prove (2).

The same statement can be phrased as a zero-energy resolvent criterion.
On the killed gap half-line, the harmonic threshold mode is `h(d)=d`.
A finite local insertion has a regular matrix element between the two
Dirichlet threshold waves. It changes their scattering length but cannot
replace either factor h by the constant bulk mode. Only a wall-bypass term
with a threshold kernel singular as

```text
[(p_2-p_1)(q_2-q_1)]^-1                                     (14)
```

could remove both external zeros and produce a nonzero
`kappa_hp(infinity)`.

## 4. Strict zero criterion

The argument proves the following reusable but narrow criterion.

> Let B be a single-insertion operator for the two directed boundaries.
> If (i) `B=P_-BP_-`, (ii) its gap/collar kernel is uniformly summable after
> stripping its declared bulk activity, and (iii) its placement is normalized
> by `1/L`, then
> `I0(2c)^(-2)<partial|D_L(Q,B)|partial>=O(c^-1)`.
> Consequently its endpoint coefficient at the unrestricted-bulk scale is
> zero.

The physical operator (8)-(10) satisfies all three conditions. Conversely,
a nonzero coefficient would identify a precise failure:

- allowing the two contours to cross during the packet;
- identifying a shared cut edge as a legal intermediate state;
- retaining an unstripped long-range gap tail; or
- summing root positions without the `1/L` already contained in beta.

None is present in the stated definition of `H_c^phys`. Therefore

```text
kappa_hp(infinity)=0                                         (15)
```

is a structural Dirichlet zero, not a numerical cancellation between
positive packet types.

## 5. Consequence for the power gate

Substituting (2) into the one-defect expansion gives

```text
beta kappa_hp(c)=O(beta/c),
[beta kappa_hp(c)]/rho(c)=O(beta).                            (16)
```

Thus the single physical hairpin cannot saturate the previously sufficient
`alpha=c beta` error bound. For this mechanism, the natural vanishing
condition is only

```text
beta=L/m^2 ->0.                                               (17)
```

The earlier `O(alpha)` estimate remains a valid relaxed absolute bound
because it forgot the antisymmetric zeros while comparing the hairpin packet
with `I0^2`. Equation (15) shows why that bound is not sharp for a correctly
normalized, physical one-hairpin insertion.

## Scientific card

- **Operator fixed:** equations (8)-(10) define the rooted physical
  single-hairpin insertion and its beta normalization.
- **Exact mechanism:** reflection keeps the evolution in the antisymmetric,
  killed-gap sector; the placement resolvent (12) is regular.
- **New value:** `kappa_hp(infinity)=0`.
- **Power consequence:** this packet changes the endpoint only by relative
  `O(beta)`, so it does not create an alpha-scale crossover.
- **Boundary:** the result covers one finite physical hairpin. It does not
  sum multiple hairpins, nonlocal collar states or an operator that actually
  opens the Dirichlet wall.
