# CM/Hecke arithmetic RG and an elliptic-point spin-4 annihilator

**Status:** C0 / deliberately speculative research proposal, 2026-08-28.

This note proposes a stronger organizing principle for the Gaussian-torus results. The claim is not that the following structure has been established. The claim is that it is sufficiently rigid to generate parameter-free or nearly parameter-free tests that are qualitatively different from another exponent fit.

## Core thesis

The Gaussian multiplier experiments may be the visible finite-lattice shadow of an **arithmetic renormalization action** at the square complex-multiplication (CM) torus.

At continuum modulus `tau=i`, the elliptic curve has endomorphism ring `Z[i]`. Multiplication by a Gaussian integer

```text
h = a + i b
Q = Norm(h) = a^2+b^2
```

is a degree-`Q` self-isogeny. On the square lattice the same operation is already used by this repository as an exact finite graph cover: it multiplies site count by `Q`, rotates microscopic orientation by `arg h`, and leaves the continuum torus modulus fixed.

The aggressive hypothesis is therefore:

> finite-size correction sectors on the square CM torus approximately furnish representations of the multiplicative semigroup `Z[i]\{0}`, and the observed Gaussian multiplier ratios are its eigenvalues.

For a semisimple scaling sector with spin `s` and radial exponent `alpha` in site count `N`, the natural complex eigencharacter is

```text
chi_(s,alpha)(h)
    = Norm(h)^(-alpha) * (h/|h|)^s.
```

For square-lattice harmonics `s in 4 Z`, the phase is insensitive to multiplication of `h` by a Gaussian unit. The angular factor is therefore the same arithmetic object that appears as the infinity-type part of a Hecke character on `Q(i)`; the extra norm power is a real quasicharacter.

The repository's real `cos(s theta)` contrasts are real projections of this complex character. The existing `1+i` and `2+-i` tests can thus be read as the first two prime-isogeny spectroscopy points, not merely two convenient geometric transformations.

This reframing becomes useful only if it predicts something new. It does, in three directions below.

---

## H1. Gaussian multiplier spectroscopy should obey exact composition

For a pure eigen-sector,

```text
chi(h1*h2) = chi(h1) chi(h2).
```

This makes unique factorization in `Z[i]` an experimental-design tool. Instead of fitting powers on unrelated sizes, build commuting diagrams whose two paths have the same endpoint.

For example,

```text
h2 = 1+i       Norm=2
h5 = 2-i       Norm=5
h10 = h2*h5 = 3+i   Norm=10
```

so a parent at `N` has a common `10N` descendant through

```text
N -> 2N -> 10N
N -> 5N -> 10N.
```

For the existing source sizes this gives

```text
65 -> 130 -> 650
65 -> 325 -> 650

85 -> 170 -> 850
85 -> 425 -> 850.
```

The endpoint is not a second measurement of the same claim; it is a **path-independence/flatness constraint** on the whole multiplier action. Any proposed mixture of H4/H12/radial sectors has to compose consistently around this square.

A useful implementation should preserve the exact Gaussian multiplier, canonicalization/unit factors, and D4 action so the complex character can be scored before taking its real projection.

### Strong version

After resolving enough independent multiplier primes, fit a small representation of the Gaussian semigroup, not separate exponents. If a putative field requires multiplier-dependent exponents or arbitrary path-dependent amplitudes, reject the representation picture even if individual pairwise fits look acceptable.

---

## H2. LCFT logarithms should appear as a non-semisimple Hecke/Jordan action

The q=2 versus Jordan-log ambiguity in `P4[S']` suggests a sharper test than another radial fit.

Suppose a rank-2 logarithmic pair sits over a leading eigencharacter `chi(h)`. Then scale multiplication should act schematically as

```text
rho(h) = chi(h) * exp(log Norm(h) * J),
J^2 = 0.
```

Equivalently, after removing the leading power/phase, a rescaled observable has

```text
Y_N = A + B log N
```

(up to the convention-dependent factor of two between `log N` and `log L`). Therefore

```text
Y_(QN) - Y_N = B log Q,
```

which is independent of the parent size. The logarithmic drift is an **additive 1-cocycle** on the multiplicative norm semigroup:

```text
delta(Q1 Q2) = delta(Q1) + delta(Q2).
```

That gives a very clean discriminator against an ordinary analytic correction. For example, if

```text
Y_N = A + C/N,
```

then the increment depends on the parent `N` and does not define an additive cocycle.

### The Q=2,5,10 plaquette curvature

Define the semigroup plaquette statistic

```text
K_N = Y_(10N) - Y_(5N) - Y_(2N) + Y_N.
```

Then, at the declared leading order,

```text
pure eigenfield:  K_N = 0
rank-2 log:       K_N = 0
A + C/N:          K_N = (2/5) C/N
```

because `1/10 - 1/5 - 1/2 + 1 = 2/5`.

This is a qualitatively different observable from selecting between `A+B/N` and `A+B log N` on a short size interval. It tests whether the correction really exponentiates as a non-semisimple representation under exact multiplier composition.

The natural first target is

```text
Y_N = N^(5/4) P4[S']
```

because the repository already has `N,2N` data and plans the `5N` child. A later `10N` common child closes the square.

A nonzero `K_N` with the sign/magnitude predicted by the frozen analytic model would be direct evidence against a pure rank-2 cocycle explanation. A near-zero `K_N` on two independent parent lineages would make the Jordan interpretation substantially more concrete than a descriptive logarithmic fit.

---

## H3. The spin-4 thermal descendant may have an E4 modular shape factor

A 2026 result changes the theory landscape enough that the torus-modulus program should be made more aggressive.

Roux, Ribault and Jacobsen, arXiv:2604.24491, derive torus one-point functions in critical loop models and show modular covariance. For the Potts energy operator `V^d_<1,2>`, fusion leaves a single torus block and they obtain

```text
<V^d_<1,2>> = |F|^2,
F = eta(q)^(2 Delta_(1,2)).
```

At percolation, `c=0` and `Delta_(1,2)=5/8`, so the chiral energy one-point block is exceptionally simple:

```text
F_energy(tau) = eta(tau)^(5/4)
```

up to the paper's nome convention.

The repository has independently identified a non-null level-4 quasiprimary in this same thermal Virasoro module,

```text
Q4 = (40 L_-2^2 - 60 L_-3 L_-1 - 9 L_-4)|h=5/8>,
```

with spin `+4` after tensoring with the antiholomorphic thermal primary (and the conjugate spin `-4` field).

Here is the strong minimal closure hypothesis:

> If the lattice matching-odd H4 residual couples to this ordinary thermal quasiprimary without an additional vector-valued/logarithmic torus sector, then the chiral ratio
>
> `F_Q4(tau) / F_energy(tau)`
>
> is a holomorphic modular form of weight 4.

The space of full-level holomorphic modular forms of weight 4 is one-dimensional. Under those assumptions,

```text
F_Q4(tau) / F_energy(tau) = C * E4(tau)
```

for a constant `C` fixed by one normalization/OPE matrix element.

This is **not** being asserted as established. Logarithmic mixing, a nontrivial defect/vector-valued modular representation, or the fact that the matching observable is not literally the Potts energy one-point function can invalidate this minimal closure. But the hypothesis has an unusually sharp consequence that does not require knowing `C`.

---

## H4. The hexagonal elliptic point is a geometric H4 annihilator

Let

```text
omega = exp(i*pi/3) = 1/2 + i*sqrt(3)/2.
```

This is `T`-equivalent to the order-3 elliptic fixed point `rho=exp(2*pi*i/3)`. The normalized Eisenstein series satisfies

```text
E4(omega) = 0,
```

with a simple zero.

The same conclusion can be viewed without naming `E4`. The stabilizer of the hexagonal torus acts locally as a `pi/3` rotation. A scalar-sector torus one-point insertion with spin increment `s` can survive the fixed point only when the stabilizer phase is trivial. Thus

```text
spin 4:   forbidden at the elliptic point
spin 12:  allowed.
```

Compare this with the square CM point `tau=i`, whose order-2 modular stabilizer corresponds to a `pi/2` rotation:

```text
tau=i:      spin 4 and spin 12 both allowed
tau=omega:  spin 4 killed, spin 12 allowed.
```

This is an **orthogonal harmonic discriminator** to norm-5. Norm-5 changes microscopic orientation at fixed square modulus. The elliptic-point experiment changes continuum modulus so that the moduli-space stabilizer itself projects out H4.

If an apparent H12 alias is merely H4 sampled at insufficient angular rank, it should disappear with the H4 sector near `omega`. If a genuine H12 component exists, the hexagonal point is precisely where it can become relatively enhanced.

---

## H5. Pell approximants give an explicit square-lattice route to the elliptic zero

The square lattice cannot realize an exactly equilateral period basis with integer vectors, but it can approach the hexagonal modulus arithmetically fast.

Take integer period vectors

```text
u = (2b, 0)
v = (b, a).
```

Their continuum modulus and site count are

```text
tau_(a,b) = 1/2 + i*a/(2b),
N = det(u,v) = 2ab.
```

Choose Pell approximants to `sqrt(3)`:

```text
a^2 - 3 b^2 = D,
```

with small fixed `D` such as `+1` or `-2`. Then

```text
tau_(a,b) - omega
  = i * (a-sqrt(3)b)/(2b)
  = i * D/[2b(a+sqrt(3)b)]
  = O(1/N).
```

More precisely,

```text
N * Im(tau_(a,b)-omega)
  = D * a/(a+sqrt(3)b)
  -> D/2.
```

Since `E4` has a simple zero,

```text
E4(tau_(a,b)) = O(1/N).
```

A direct q-series evaluation gives the useful geometric diagnostics below. These numbers are **modular-shape proxies, not fitted percolation amplitudes**.

| `(a,b)` | `D=a^2-3b^2` | `N=2ab` | `Im tau` | `E4(tau)` | `N E4(tau)` |
|---|---:|---:|---:|---:|---:|
| `(5,3)` | -2 | 30 | 0.8333333333 | -0.2170053400 | -6.51016 |
| `(7,4)` | +1 | 56 | 0.8750000000 | +0.0527808750 | +2.95573 |
| `(19,11)` | -2 | 418 | 0.8636363636 | -0.0145179092 | -6.06849 |
| `(26,15)` | +1 | 780 | 0.8666666667 | +0.00386292016 | +3.01308 |
| `(71,41)` | -2 | 5822 | 0.8658536585 | -0.00103701130 | -6.03748 |
| `(97,56)` | +1 | 10864 | 0.8660714286 | +0.000277726626 | +3.01722 |

The `D=+1` and `D=-2` families approach opposite sides of the zero, with asymptotic proxy amplitudes in an approximately `-2:1` ratio as expected from the Pell defect `D`.

This gives two built-in controls: convergence rate and sign reversal are fixed by geometry before percolation data are generated.

---

## H6. Geometric superconvergence: L^-4 root bias -> L^-6

The current leading matching-odd H4 picture is

```text
M_H4(pc; tau=i) ~ N^(-13/8)
M'(pc)          ~ N^(+3/8),
```

which gives the familiar root bias

```text
Delta p ~ N^-2 = L^-4.
```

Under the minimal modular closure, the H4 amplitude near the elliptic point acquires the simple-zero factor

```text
E4(tau_N) ~ N^-1.
```

Along a fixed-`D` Pell family,

```text
M_H4(pc; tau_N) ~ N^(-13/8) * N^-1
                = N^(-21/8).
```

Dividing by the thermal slope predicts

```text
Delta p_H4 ~ N^(-21/8-3/8)
           = N^-3
           = L^-6.
```

So the bold prediction is:

> A Pell-to-hexagonal torus sequence should geometrically improve the H4-induced pseudo-critical root convergence from `L^-4` to `L^-6`, without a multi-size annihilator, if H4 is genuinely controlled by a simple modular zero and no lower surviving sector takes over.

This is potentially more interesting than another threshold estimate. It would be a finite-size **improved geometry** produced by a zero of a continuum shape amplitude.

### What failure would mean

Failure is highly diagnostic rather than generic:

- residual remains `N^-13/8`: a surviving H12/higher harmonic or non-H4 sector is dominating;
- H4 is suppressed but not as `1/N`: the shape amplitude is not the minimal `E4` closure or the integer-torus approach introduces another scaling field;
- root improves but to a different power: another correction becomes leading after H4 annihilation;
- opposite Pell families do not reverse the leading signed shape residual: the simple-zero modular picture is wrong or the measured observable is in a vector-valued/logarithmic sector.

---

## H7. Two CM points may form a useful operator spectroscope

The square and hexagonal elliptic points should be treated as complementary high-symmetry laboratories:

```text
Gaussian CM point:   tau=i,      End ~ Z[i]
Eisenstein CM point: tau=omega,  End ~ Z[omega]
```

At `tau=i`, Gaussian integer self-isogenies provide exact scale+rotation eigenvalue tests. Near `tau=omega`, the larger stabilizer supplies spin selection rules and an H4 zero.

This suggests a broader program:

> finite-size correction amplitudes are vector-valued modular objects, while exact lattice covers sample arithmetic correspondences on moduli space; semisimple irrelevant fields appear as multiplier eigencharacters and LCFT partners as non-semisimple/Jordan extensions of those characters.

In that language, the present `Delta cos(4 theta) N^-13/8` law is the first observed matrix element of an arithmetic RG representation.

This may be wrong. But it is a much stronger hypothesis than “there is a power law”, and it organizes angular spin, radial exponent, logarithmic mixing, torus modulus and exact Gaussian covers in one algebraic object.

---

## Concrete research order

This proposal should not displace already-running target data. It changes what to do with the next expensive geometry.

### 1. Zero-extra-compute arithmetic reanalysis

For every existing Gaussian multiplier result, store the full complex multiplier `h` and test multiplicative character consistency rather than only norm/angular scalar ratios.

For P48-style rescaled channels, report semigroup increments

```text
Y_(QN)-Y_N
```

for every available `Q` and ask whether they are closer to parent-independent `B log Q` or to analytic parent-dependent corrections.

### 2. Close one norm-10 commuting square after norm-5

If N=325/425 production is scientifically successful, do not automatically buy another unrelated large size. Estimate the information-per-CPU of `N=650` and/or `N=850` because each closes a `2 x 5` multiplier plaquette and directly tests representation composition/Jordan curvature.

### 3. Implement the cheap end of the Pell elliptic sequence

The general integer-period homology engine should be able to represent

```text
u=(2b,0), v=(b,a).
```

Start with exact/small controls and then the moderate pair

```text
D=+1: (a,b)=(7,4),   N=56
D=-2: (a,b)=(19,11), N=418
```

or choose a better cost-matched pair after variance profiling. The important design feature is **opposite sides of the same modular zero**, not equal N.

### 4. Freeze an elliptic-zero score before large runs

A useful primary score is not an absolute amplitude. Use a ratio/normalized quantity in which the unknown microscopic H4 coupling cancels as much as possible. At minimum score:

```text
N * A4_shape(tau_N)
```

for stability within a fixed Pell `D` family and sign reversal between `D=+1` and `D=-2`.

If a derivation of the descendant torus one-point ratio fixes the proportionality to `E4`, promote the full `E4(tau)` shape ratio to a parameter-free target.

### 5. Search what is revealed after H4 is killed

If the elliptic annihilator works, immediately inspect the residual for:

- H12/higher angular content;
- scalar matching-odd sectors;
- the post-H4 correction exponent behind the historical annihilator;
- logarithmic shape response.

An improved geometry is valuable because it turns a dominant correction into a controlled zero and exposes the next field.

---

## Theory calculations worth doing before more model proliferation

1. **Descendant one-point calculation.** Starting from the one-block Potts energy torus amplitude, use Zhu/Virasoro recursion to compute the torus one-point function of the explicit `Q4` state. Determine whether the ratio to `eta^(5/4)` is exactly proportional to `E4`, or whether a logarithmic/vector-valued term is forced.
2. **Matching-observable bridge.** Determine whether the FK/topological-sector combination representing the matching function couples to the same one-dimensional energy block or to a larger modular representation.
3. **Jordan multiplier algebra.** Write the minimal 2x2 semigroup representation and derive joint Q=2/Q=5/Q=10 covariance-aware residuals for `P4[S']`.
4. **Elliptic-point selection table.** For every currently live H0/H4/H8/H12 candidate, derive the stabilizer selection rule at `tau=i` and `tau=omega`, including possible vector-valued phases.
5. **CM-to-CM extension.** Ask whether an exactly critical control model on an isoradial/triangular representation can realize the Eisenstein CM point exactly and calibrate the predicted spin annihilation before applying it to square-site approximants.

---

## Relation to existing issues

This is related to, but stronger/more specific than:

- #57 / #64: Gaussian multiplier spectroscopy — reinterpreted as semigroup/Hecke eigenvalue spectroscopy and extended to composition curvature;
- #103: torus-modulus spectroscopy — upgraded with a concrete elliptic-point zero and Pell sequence;
- #125: operator mixing — the proposed semigroup representation supplies a constrained mixing algebra;
- #48/P48: q=2 versus Jordan drift — replaced by a cocycle/plaquette test that directly targets non-semisimplicity;
- #37: `x=21/4` thermal spin-4 candidate — supplied with a new modular-shape fingerprint;
- #106: improved-action controls — here the improvement is geometric/modular rather than a tuned microscopic coupling.

The proposal is intentionally not a new claim in `STATUS.md`. It is a research program that becomes valuable precisely because several of its predictions can fail cleanly.

---

## External anchors

- P. Roux, S. Ribault, J. L. Jacobsen, **Torus one-point functions in critical loop models**, arXiv:2604.24491 (2026). In particular Sec. 5.1 gives the one-block Potts energy torus one-point function and modular covariance.
- M. R. Gaberdiel, S. Lang, **Modular differential equations for torus one-point functions**, arXiv:0810.0106.
- A. Poghosyan, **Shaping Lattice through irrelevant perturbation: Ising model**, arXiv:1908.06291. This is an explicit precedent for identifying square-lattice finite-size corrections with irrelevant CFT perturbations on the torus.
- L.-P. Arguin, **Homology of Fortuin-Kasteleyn clusters of Potts models on the torus**, arXiv:hep-th/0111193.
- J. Dubail, J. L. Jacobsen, H. Saleur and related periodic Potts/Temperley-Lieb literature; see also the repository's existing torus-sector references.
- Standard CM fact: the square elliptic curve has endomorphism ring `Z[i]`, while the hexagonal elliptic curve has endomorphism ring `Z[omega]`; endomorphism degree is the algebraic norm.
- Standard modular-form facts: `M_4(SL(2,Z))` is one-dimensional and generated by `E4`; `E4` has its elliptic zero at the order-3 point, while `E6` vanishes at the order-2 point.

## Bottom line

The strongest speculative synthesis is now:

```text
Gaussian covers
  -> CM self-isogenies at tau=i
  -> multiplicative scaling/spin characters
  -> non-semisimple Jordan extensions for LCFT logs
  -> commuting-square curvature tests

modulus deformation
  -> elliptic stabilizer spin selection
  -> E4 zero at tau=exp(i*pi/3)
  -> Pell O(1/N) approach from integer square-lattice tori
  -> H4 residual N^(-13/8) -> N^(-21/8)
  -> root bias N^-2 -> N^-3 (L^-4 -> L^-6)
```

If even half of this survives, the project stops being only a finite-size percolation study and becomes an arithmetic/modular spectroscopy program for lattice corrections at `c=0`.