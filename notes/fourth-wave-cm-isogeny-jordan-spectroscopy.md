# Fourth-wave research: Gaussian CM spectroscopy, Jordan cocycles, and improved percolation

Status: **C0 research program / deliberately aggressive theory note**.

Related issue: #138.

This note starts from the current numerical picture rather than re-arguing it. The repository already has a reproducible matching-odd square-harmonic sector, approximate `N^-13/8` radial behavior, three successful `1+i` Gaussian lineages, a resolved local residual-to-root mechanism, and a derivative channel `P4[S']` whose pure `N^-5/4` law is insufficient at current precision.

The next conceptual step is to stop treating Gaussian multiplication as a collection of clever finite-size ratios and instead ask:

> **What representation of the Gaussian cover semigroup is carried by the finite-size correction sector?**

For an ordinary CFT eigenfield the answer should be a multiplicative character. For an LCFT Jordan block it should be a generalized character with logarithmic cocycles. This viewpoint unifies radial exponent, conformal spin, Gaussian genealogy, logarithmic mixing, and future experiment design.

The strongest new practical consequence is immediate: the existing norm-2 full curves plus the already-planned norm-5 children give a parameter-free `q=2` versus Jordan discriminator without fitting the correction amplitudes at all.

---

## 1. The square Gaussian torus is a CM object, not just a geometry trick

The continuum square torus is

\[
E_i=\mathbb C/(\mathbb Z+i\mathbb Z).
\]

Its endomorphism ring contains the Gaussian integers. For

\[
h=a+ib\in\mathbb Z[i],
\]

multiplication

\[
z\mapsto hz
\]

is a complex-multiplication endomorphism of degree

\[
Q=N(h)=a^2+b^2=|h|^2.
\]

On the repository's primitive Gaussian quotients the same arithmetic operation is visible discretely: it multiplies the site count by `Q`, physical length by `|h|`, and microscopic orientation by

\[
\phi_h=\arg h,
\]

while preserving the square modular shape.

So the current exact graph-cover hierarchy is the finite-lattice shadow of a CM endomorphism hierarchy of the square torus.

This changes how I would organize the experiment space. `Norm=2`, `Norm=5`, `Norm=13`, ... should not be thought of as unrelated scale factors. They are arithmetic probes of one multiplicative action.

### Prime arithmetic gives a natural spectroscopy ladder

The first cases are structurally special:

- `2=(1+i)(1-i)` is ramified in `Z[i]`. The angle is `pi/4`, which aliases all odd square harmonics `H4,H12,H20,...` by the same sign reversal.
- `5=(2+i)(2-i)` is the first split odd prime. It supplies a nontrivial orientation-changing degree-5 cover and breaks the norm-2 H4/H12 alias.
- more generally, primes `p=1 mod 4` split and furnish non-axis degree-`p` Gaussian covers;
- primes `p=3 mod 4` are inert and do not furnish a non-axis Gaussian element of norm `p`.

Thus #57 is not merely "the cheapest multiplier after 2". It is the **first split-prime CM spectroscopy experiment**.

This suggests a future design rule: after norm 5, choose split primes by angular-model information gain per CPU cost, not by increasing `N` monotonically.

---

## 2. Pure correction fields should carry Gaussian multiplicative characters

Consider a dimensionless torus observable with a correction from a continuum field of scaling dimension `x` and spin `s`. In a fixed physical frame write the complex harmonic form

\[
\delta X_{x,s}(N,\theta)
\sim
N^{-\alpha}\,\operatorname{Re}\!\left(A_s e^{is\theta}\right),
\qquad
\alpha=\frac{x-2}{2}.
\]

Under multiplication by `h`,

\[
N\mapsto QN,
\qquad
\theta\mapsto\theta+\phi_h.
\]

Therefore a complex spin projection should obey

\[
\boxed{
\mathcal P_s(QN)
=
Q^{-\alpha}e^{is\phi_h}\mathcal P_s(N)
+	ext{subleading sectors}.
}
\]

Define

\[
\lambda_{x,s}(h)
=|h|^{2-x}e^{is\arg h}
=Q^{-\alpha}e^{is\phi_h}.
\]

Then at the pure-field level

\[
\boxed{
\lambda_{x,s}(h_1h_2)
=\lambda_{x,s}(h_1)\lambda_{x,s}(h_2).
}
\]

That is the core conjecture: **a correction eigenfield furnishes a one-dimensional representation of the Gaussian multiplicative semigroup.**

The current no-fit doubling laws are exactly what one expects from this structure after taking a real two-angle projection.

### Current candidate sectors in this language

For the matching-even square anisotropy candidate

\[
x_I=4,\quad s=4,
\]

we have

\[
\alpha_I=1.
\]

For the matching-odd thermal-family candidate

\[
x_T=21/4,\quad s=4,
\]

we have

\[
\alpha_T=13/8.
\]

The first thermal derivative multiplies a projected finite-size term by `N^(3/8)`, producing the repository's `5/4` exponent in `P4[S']`.

### Why the current real contrast is only a shadow of the representation

For a two-angle real contrast, Gaussian multiplication gives

\[
\frac{\Delta X_{QN}}{\Delta X_N}
=
Q^{-\alpha}
\frac{\Delta\cos[s(\theta+\phi_h)]}{\Delta\cos(s\theta)}
\]

when the lattice coupling is reflection-even in the chosen frame. The pair-dependent angular ratio is why the exact norm-5 raw H4 transfer is not simply `cos(4 phi_h)`.

A more fundamental object would reconstruct both cosine and sine components and let the multiplier act by the literal phase `exp(i s phi_h)`. The N=1105 multi-angle machinery is already close to what is needed for this complex spin tomography.

---

## 3. The LCFT upgrade: Jordan fields are generalized Gaussian eigenvectors

Percolation is logarithmic. Vasseur--Jacobsen--Saleur identify logarithmic mixing involving the energy operator at `Q=1`, and recent work by Yifei He makes the stronger structural point that Kac operators at `c=0` sit at the bottom of logarithmic multiplets and that higher-rank Jordan structures are plausible.

This suggests replacing the scalar character above by a triangular representation.

For a rank-2 Jordan pair, scaling by a factor `ell` acts schematically as

\[
\begin{pmatrix}
\phi\\
\psi
\end{pmatrix}
\mapsto
\ell^{2-x}
\begin{pmatrix}
1&0\\
\gamma\log\ell&1
\end{pmatrix}
\begin{pmatrix}
\phi\\
\psi
\end{pmatrix}.
\]

Adding spin under a Gaussian multiplier gives

\[
\rho(h)
=
|h|^{2-x}e^{is\arg h}
\exp\!\left[\Gamma\log|h|\right],
\]

where `Gamma` is nilpotent inside the Jordan block.

The important point is not notation. It is that Gaussian multiplication still composes exactly:

\[
\rho(h_1h_2)=\rho(h_1)\rho(h_2),
\]

but the representation is no longer diagonalizable.

A logarithm is therefore not merely an optional correction term in a fit. It is the observable cocycle of a generalized eigenvector under the same semigroup already being probed by norm-2 and norm-5 experiments.

---

## 4. A new zero-amplitude-fit discriminator using N, 2N, and 5N

This is the most immediately useful consequence of the framework.

For a normalized projected channel with frozen leading exponent `alpha`, define

\[
Y(N)=N^\alpha P_s[X](N).
\]

The live `P4[S']` mechanisms have qualitatively different dependence on multiplicative scale.

### 4.1 Rank-2 Jordan law

If

\[
Y(N)=A+B\log N,
\]

then

\[
Y(QN)-Y(N)=B\log Q.
\]

For the same parent with norm-2 and norm-5 descendants,

\[
\boxed{
\mathcal C_{5|2}(N)
=
\frac{Y(5N)-Y(N)}{Y(2N)-Y(N)}
=
\frac{\log5}{\log2}
=2.321928094887362\ldots
}
\]

The amplitude `A`, logarithmic coefficient `B`, and parent size all cancel.

### 4.2 Ordinary relative-power law

If instead

\[
Y(N)=A+C N^{-\beta},
\qquad \beta=q/2,
\]

where `q` is the correction exponent in physical length, then

\[
\boxed{
\mathcal C_{5|2}(q)
=
\frac{5^{-q/2}-1}{2^{-q/2}-1}.
}
\]

Fixed targets are

```text
relative length q      C_5|2
--------------------------------
1                      1.88733083934
3/2                    1.72899949352
2                      1.60000000000
3                      1.40855759416
4                      1.28000000000
6                      1.13371428571
rank-2 Jordan log       2.32192809489
```

The `q=2` target is especially important because it is the already-frozen ordinary-RG correction competing with the Jordan model in P48.

The `q=3/2` row is worth retaining as a generic-percolation adversary: Xu--Chen--Zhou--Salas--Deng derived the cluster-size correction exponent `Omega=72/91`, which converts through the percolation fractal dimension `D=91/48` to a length correction `omega=D Omega=3/2`. This does **not** prove that the same field appears in the projected torus derivative channel; it gives a modern theory-motivated leakage candidate rather than an arbitrary fractional power.

### 4.3 Use a linear residual, not a noisy literal ratio

For a frozen candidate `c`, score

\[
\boxed{
R_c(N)
=
[Y(5N)-Y(N)]
-c[Y(2N)-Y(N)].
}
\]

Equivalently,

\[
R_c=Y(5N)-cY(2N)+(c-1)Y(N).
\]

This is linear in the measured points, so covariance propagation is straightforward and the statistic remains well-defined if one increment is small.

### 4.4 The experiment is already almost paid for

For `P4[S']`, use

\[
\alpha=5/4.
\]

The two triples are

```text
N=65:  65 -> 130 via norm 2; 65 -> 325 via norm 5
N=85:  85 -> 170 via norm 2; 85 -> 425 via norm 5
```

The norm-2 full curves already exist. Issue #57 is already planning the N=325 and N=425 full curves.

So after #57 there is no reason to compare `q=2` and Jordan only through separately fitted two-parameter curves. We can ask the sharper question:

> Does the measured scale increment behave like a power character or a logarithmic cocycle?

The corresponding frozen artifact is `predictions/p48_sprime_semigroup_curvature_20260828.yaml`.

---

## 5. A norm-10 commuting square would measure Jordan rank, not just another exponent

If the three-point curvature test remains ambiguous, the next child should not be chosen merely because it is larger. Use the product of the two existing multipliers.

For example,

\[
(1+i)(2-i)=3+i,
\qquad N(3+i)=10.
\]

This yields a multiplicative square

```text
       norm 2
   N ----------> 2N
   |              |
 n5|              |n5
   v              v
  5N ----------> 10N
       norm 2
```

The two paths are the same Gaussian multiplication because `Z[i]` is commutative.

For the scaled observable define the mixed multiplicative difference

\[
\boxed{
\mathcal H_{2,5}Y(N)
=Y(10N)-Y(5N)-Y(2N)+Y(N).
}
\]

Then:

### Pure eigenfield

\[
Y=A
\quad\Longrightarrow\quad
\mathcal H_{2,5}Y=0.
\]

### Rank-2 Jordan block

\[
Y=A+B\log N
\quad\Longrightarrow\quad
\mathcal H_{2,5}Y=0.
\]

### Ordinary power correction

\[
Y=A+C N^{-\beta}
\]

gives

\[
\boxed{
\mathcal H_{2,5}Y
=C N^{-\beta}(1-2^{-\beta})(1-5^{-\beta}).
}
\]

### Rank-3 Jordan block

If

\[
Y=A+B\log N+C(\log N)^2,
\]

then

\[
\boxed{
\mathcal H_{2,5}Y
=2C\log2\log5.
}
\]

The quadratic logarithmic coefficient is isolated as a scale-independent mixed curvature.

This turns potential N=650/N=850 runs into **Jordan-rank tomography**. Given the recent c=0 literature, that is a much better scientific reason to run them than "extend the radial range".

The exact graph-cover machinery in #67 is naturally compatible with this square and may eventually make a multilevel coupled estimator possible on all four corners.

---

## 6. Matching parity may be inherited by an entire Virasoro module

The repository is correctly cautious about promoting the exact lattice matching relation to a full local OPE automorphism. I think there is a stronger intermediate statement worth attacking.

Treat the square-site lattice and its matching partner as two UV regularizations of the same continuum fixed point. Write schematically

\[
S_G
=S_*+t\int\epsilon+\sum_a u_a\int O_a,
\]

\[
S_{\hat G}
=S_*-t\int\epsilon+\sum_a \hat u_a\int O_a.
\]

Occupation complement reverses the thermal direction, so the thermal primary is naturally exchange-odd in the doubled UV theory space.

Now make the **module-intertwiner hypothesis**:

> on a nondegenerate Virasoro or generalized Virasoro block, matching-pair exchange commutes with the geometric conformal generators.

If so,

\[
\mathcal T L_{-n}=L_{-n}\mathcal T
\]

inside that block. Therefore if

\[
\mathcal T|\epsilon\rangle=-|\epsilon\rangle,
\]

then every ordinary descendant satisfies

\[
\mathcal T\,U(L_{-n})|\epsilon\rangle
=-U(L_{-n})|\epsilon\rangle.
\]

In particular, the surviving level-4 spin-4 thermal quasiprimary inherits matching-odd parity.

This is much weaker than asserting that matching preserves every OPE coefficient. It only asks for an intertwiner on the scaling module. It is also much stronger than an empirical sign label because it predicts parity for the full descendant tower.

### Logarithmic refinement

In a degenerate/Jordan block, `T` need not be a scalar sign on a chosen basis. It can be a finite matrix that commutes with the nilpotent dilation operator. The correct object may therefore be a **matching-parity generalized block**, not an individually signed field.

This is another reason to elevate the Gaussian generalized-eigenvector tests: they can reveal a matrix-valued scaling block before a full OPE construction exists.

---

## 7. The April-2026 torus-one-point work opens a direct operator-fingerprint route

A particularly important literature development is Roux--Ribault--Jacobsen, *Torus one-point functions in critical loop models* (`arXiv:2604.24491`). They systematically compute torus one-point functions in critical loop models, including Potts-type models, and relate modular covariance on the torus to sphere four-point bootstrap data. Their conformal-block treatment includes logarithmic blocks.

This changes the long-term theory target for Matching One.

The current `x=21/4,s=4` claim is based mainly on:

1. radial dimension matching;
2. spin/H4 behavior;
3. an explicit level-4 Virasoro state surviving the thermal null quotient;
4. matching/thermal derivative parity.

A much stronger identification would predict the **function of torus modulus**.

### Proposed analytic pipeline

Let

\[
F_\epsilon(\tau,\bar\tau)
=\langle\epsilon\rangle_{T^2}
\]

or the appropriate homology/combinatorial-map-resolved loop-model one-point object.

The repository already has an explicit chiral level-4 quasiprimary

\[
Q_4
=(40L_{-2}^2-60L_{-3}L_{-1}-9L_{-4})|h=5/8\rangle.
\]

General torus Ward/Zhu recursion expresses Virasoro-descendant one-point functions as quasi-modular differential operators acting on primary one-point functions. Therefore the target is conceptually

\[
F_{Q_4}(\tau)
=\mathcal D_4[\tau;E_2,E_4,\ldots]F_\epsilon(\tau),
\]

with the logarithmic completion treated in the `c=0` theory.

Then one normalization at `tau=i` would predict amplitude ratios at other sheared tori:

\[
\frac{A_T(\tau_1)}{A_T(i)}
=\frac{F_{Q_4}(\tau_1)}{F_{Q_4}(i)}
\]

up to the precise frame/spin convention.

That would be an operator fingerprint. H4, H12, a different Kac family, and a generic lattice artifact can share a radial exponent over a short size range; they should not generically share the same full modular-shape function.

### Concrete research task

Rather than starting from a formal LCFT construction from scratch:

1. inspect the 2026 loop-model torus one-point solutions and code;
2. locate the Potts/percolation specialization and the field closest to the thermal Kac primary;
3. analytically continue or take the `Q->1` limit with the correct normalization;
4. apply level-4 descendant recursion;
5. generate a numerical table of the predicted shape ratio along a one-parameter torus path;
6. compare against square-site, square-bond, and self-matching controls.

This is, in my view, a higher-value theory project than another independent derivation of the exponent `21/4`.

---

## 8. Upgrade the magic-torus idea from a zero hunt to a modular-function test

The existing `magic-torus-modular-conjecture.md` correctly treats the `E4(tau_hex)=0` observation as heuristic. The new one-point-function perspective suggests a way to sharpen it.

For the matching-even identity-family spin-4 sector, the relevant stress-tensor/KdV descendants have torus expectation values controlled by modular or quasi-modular forms. In ordinary CFT, thermal KdV correlators are quasi-modular differential operators on the torus partition function.

So instead of asking only

> is there a small amplitude near the hexagonal point?

ask

> after fixing the correct quasi-modular basis, does the whole measured shape dependence lie in the finite-dimensional function space predicted for the candidate descendant?

The strongest simple model is still that the leading holomorphic weight-4 piece is proportional to `E4(tau)`. If that component dominates, the equianharmonic zero is automatic. But the decisive test should fit/predict the entire path in `tau`, including signs and ratios on both sides of the zero, not merely one near-cancellation.

For the thermal-family spin-4 candidate, apply the same philosophy to the descendant of the thermal one-point function rather than borrowing the identity-family `E4` ansatz.

---

## 9. Build an "improved percolation" model that tunes the H4 coupling to zero

There is another way to learn the operator content: change the microscopic action deliberately.

In lattice field theory and critical lattice models, one often tunes microscopic couplings so the leading irrelevant scaling field has zero amplitude. The resulting improved model exposes subleading operators much sooner.

Matching One should do the same.

### 9.1 First control: C4 FK/bond family with axis and diagonal couplings

Consider a C4-symmetric short-range family with separate axis and diagonal bond weights/couplings. The continuum small-momentum expansion has a fourth angular moment of the schematic form

\[
g_4^{\rm bare}
\propto
\sum_e w_e |e|^4\cos(4\phi_e).
\]

Axis and diagonal directions contribute opposite signs because

\[
\cos(4\cdot0)=+1,
\qquad
\cos(4\cdot\pi/4)=-1.
\]

So there should generically be a tuning along the critical manifold where the leading square H4 artifact crosses zero.

The naive fourth-moment estimate gives an axis/diagonal weight ratio of order `4:1` because diagonal bonds have `|e|^4=4`; this is only a tree-level locator, not a proposed exact critical ratio.

### 9.2 Why zero crossings are operator spectroscopy

Suppose a one-parameter microscopic family has two measured spin-4 amplitudes

\[
A_I(\lambda)
\]

for the matching-even/identity-like sector and

\[
A_T(\lambda)
\]

for the matching-odd/thermal-like sector.

If they correspond to genuinely different continuum operators, their lattice couplings are independent analytic functions of `lambda`. Generically their zeros need not coincide:

\[
\lambda_I^*\ne\lambda_T^*.
\]

A resolved separation of the two zero crossings would be direct evidence that the observed `N^-1` and `N^-13/8` H4 sectors are distinct scaling fields, not one anisotropy term appearing differently in `S` and `D`.

Conversely, a robust common zero would reveal an unexpected microscopic relation between the two couplings and would itself demand explanation.

This is an experimental version of operator mixing theory: tune the UV action and watch continuum amplitudes rotate.

### 9.3 Later target: improved matching-pair family

The bond/FK family is the easiest positive control. The more ambitious construction is a continuously tunable short-range site/hyperedge family that comes with a defined matching partner and preserves the exact pair exchange. Such a family would let us trace the matching-odd thermal spin-4 coupling itself through zero.

That could become a much cleaner route to subleading matching-odd sectors than heroic increases in `N` on the original square-site action.

---

## 10. Preserve full homology: build complex spin tomography, not only coarse wrapping channels

The repository has learned that several coarse wrapping-difference channels are configuration-identical. That does not mean topology is exhausted; it means the coarse projection has collapsed information.

Torus FK/percolation theory naturally resolves primitive winding classes

\[
\{a,b\}.
\]

Morin-Duchesne and Saint-Aubin explicitly studied the modular dependence and critical exponents of these homology probabilities, and modern Potts/loop work uses twisted/topological sectors extensively.

I would therefore preserve or reconstruct the full primitive winding label whenever computationally feasible.

The new target is a complex spin-sensitive topological moment rather than another scalar wrapping probability. Schematically one can build combinations whose transformation under a rotation/shear resolves both cosine and sine components of a spin-`s` response.

With enough independent winding/orientation components, the Gaussian multiplier should act by the full complex phase

\[
e^{is\phi_h},
\]

not merely by a sign or a real angular ratio.

This would make H4/H12 mixture estimation a linear representation problem and would provide phase-sensitive checks of Gaussian composition.

---

## 11. The generic `omega=3/2` correction is now a useful selection-rule probe

A 2025 result derives the percolation cluster-size correction exponent

\[
\Omega=72/91.
\]

Using the cluster fractal dimension

\[
D=91/48,
\]

the corresponding length exponent is

\[
\omega=D\Omega=3/2.
\]

This does not mean every torus observable must contain an `L^-3/2` term. Quite the opposite: the current projectors may suppress it by topology, matching parity, spin, or observable selection rules.

That gives a sharper interpretation of absence:

> if generic unprojected quantities show the `3/2` sector while the matching-odd H4 projector does not, the projector is behaving like an automatically improved observable.

The experiment should therefore compare the same full-curve data in:

- scalar/unprojected channels;
- matching-even H4;
- matching-odd H4;
- derivative channels.

A sector-by-sector presence/absence table for the `3/2` correction would be more informative than adding it as one more free term everywhere.

---

## 12. What I would do next, in order

### P0 — score multiplicative curvature immediately after #57

No new production beyond the already-planned N=325/425 full curves.

For both parent triples, reconstruct

```text
Y(N) = N^(5/4) P4[S'](N)
```

and score the linear residuals for

```text
q=2        c = 1.6
Jordan     c = log(5)/log(2)
q=3/2      c = 1.72899949352   [secondary generic-percolation adversary]
```

plus other already-declared correction powers only if chronologically appropriate.

This should become part of #57/#64 scoring rather than a new expensive campaign.

### P0 theory — derive the Gaussian generalized-character algebra cleanly

Write the exact transformation of:

- complex H4/H8/H12 projectors;
- matching-even/odd blocks;
- thermal derivatives;
- rank-2 and rank-3 Jordan blocks;
- canonical D4 reflections/units used by the finite Gaussian quotient implementation.

The deliverable should be an exact symbolic table for arbitrary `h in Z[i]` and a design criterion based on model separation.

### P0/P1 theory — attack the 2026 torus-one-point bridge

Try to turn the level-4 thermal candidate into a predicted shape function. This is the most promising path from "candidate dimension" to "operator fingerprint" that I found in the current literature.

### P1 — only if needed, run the norm-10 commuting square

Use N=650 and N=850 (or a cheaper equivalent if the design optimizer finds one) to distinguish ordinary powers from Jordan rank and to test Gaussian composition closure.

### P1 — implement an improved square-lattice control

Tune a C4 axis/diagonal FK model through a spin-4 zero and measure how the known identity-like anisotropy disappears. Use this to validate the amplitude-zero methodology before attempting a matching-pair improved family.

### P1/P2 — full homology complex-spin archive

Before future production engines throw away winding detail, estimate the storage/runtime cost of retaining primitive homology labels or sufficient Fourier moments of them.

---

## 13. Strong conjectures worth trying to kill

I would explicitly keep the following on the board.

### Conjecture G1 — Gaussian character law

The leading projected finite-size fields organize into finite-dimensional representations of the Gaussian cover semigroup, with ordinary fields diagonal and LCFT blocks upper-triangular.

### Conjecture G2 — thermal-module parity inheritance

Matching pair exchange acts as an intertwiner on the thermal generalized Virasoro module; the level-4 spin-4 descendant therefore inherits the thermal primary's matching-odd block.

### Conjecture G3 — `P4[S']` is a Jordan cocycle

The observed derivative drift is dominated by a rank-2 logarithmic partner of the same `x=21/4,s=4` thermal sector, and the N/2N/5N curvature tends to `log5/log2` rather than `1.6`.

### Conjecture G4 — a higher-rank logarithmic sector is experimentally visible

If N/2N/5N rejects both a single power and a rank-2 log, a norm-10 mixed difference will reveal a nonzero scale-independent component consistent with a `(log N)^2` term.

### Conjecture G5 — the thermal spin-4 torus shape is bootstrap-computable

The 2026 critical-loop torus one-point framework plus descendant recursion can predict the `tau` dependence of the candidate `x=21/4` field up to one microscopic coupling.

### Conjecture G6 — there exists an improved square action with vanishing leading H4 coupling

A short-range C4-symmetric tunable model has a critical point/manifold where the leading identity-family spin-4 lattice artifact is zero, and an analogous matching-pair family can eventually tune the thermal H4 coupling independently.

### Conjecture G7 — the matching projector suppresses the generic `omega=3/2` sector

The recently derived generic percolation correction is visible in less-projected observables but absent or strongly suppressed in the matching-odd H4 sector because of combined topology/parity/spin selection.

---

## 14. Literature anchors and what they contribute

- X. Feng, Y. Deng, H. W. J. Blöte, *Percolation transitions in two dimensions*, Phys. Rev. E 78, 031136 (2008), arXiv:0901.1370. Establishes the importance of `X_t2=4`, logarithmic companions, and strong orientation dependence of the power-law correction.
- R. Vasseur, J. L. Jacobsen, H. Saleur, *Logarithmic observables in critical percolation*, J. Stat. Mech. L07001 (2012), arXiv:1206.2312. Gives the concrete energy/log-partner collision at `Q=1`.
- Y. He, *Logarithmic operators in c=0 bulk CFTs*, SciPost Phys. 19, 008 (2025), arXiv:2411.18696. Places Kac operators at the bottom of logarithmic multiplets and motivates higher-rank Jordan possibilities.
- M. Caselle, M. Hasenbusch, A. Pelissetto, E. Vicari, *Irrelevant operators in the two-dimensional Ising model*, J. Phys. A 35, 4861 (2002), arXiv:cond-mat/0106372. A useful solved-model analogue for square-lattice nonzero-spin irrelevant operators and amplitude cancellations.
- Y. Xu, T. Chen, Z. Zhou, J. Salas, Y. Deng, *Correction-to-scaling exponent for percolation and the Fortuin--Kasteleyn Potts model in two dimensions*, Phys. Rev. E 111, 034108 (2025), arXiv:2411.12646. Derives `Omega=72/91`, corresponding to a length correction `omega=3/2` for percolation clusters.
- R. M. Ziff, C. D. Lorenz, P. Kleban, *Shape-dependent universality in percolation*, Physica A 266, 17 (1999), arXiv:cond-mat/9811122. Establishes torus aspect-ratio/twist dependence as a universal continuum datum.
- A. Morin-Duchesne, Y. Saint-Aubin, *Critical exponents for the homology of Fortuin--Kasteleyn clusters on a torus*, Phys. Rev. E 80, 021130 (2009), arXiv:0812.2925. Provides a natural full-homology/twist language beyond coarse wrap events.
- M. R. Gaberdiel, S. Lang, *Modular differential equations for torus one-point functions*, J. Phys. A 42, 045405 (2009), arXiv:0810.0106. Demonstrates the modular-differential control of torus one-point functions in Virasoro theories.
- A. Maloney, G. S. Ng, S. F. Ross, I. Tsiares, *Thermal correlation functions of KdV charges in 2D CFT*, JHEP 02 (2019) 044, arXiv:1810.11053. Shows that stress-tensor/KdV torus observables are governed by quasi-modular differential operators.
- P. Roux, S. Ribault, J. L. Jacobsen, *Torus one-point functions in critical loop models*, arXiv:2604.24491 (2026). The most directly relevant new analytic development: systematic loop/Potts torus one-point functions, modular covariance, and logarithmic conformal blocks.
- J. A. Harvey, Y. Wu, *Hecke relations in rational conformal field theory*, JHEP 09 (2018) 032, arXiv:1804.06860. Conceptual precedent for Hecke/isogeny structures acting on modular CFT data. I am **not** claiming the present observable is a Hecke eigenform; the useful analogy is that isogeny composition can organize spectral data.

For the arithmetic geometry itself, the key elementary fact is that the square elliptic curve/complex torus has CM by `Z[i]`, and multiplication by a Gaussian integer is an isogeny whose degree is its norm. Hecke correspondences, by contrast, average over index-`n` sublattices/degree-`n` isogenies. The repository currently follows selected CM endomorphisms rather than performing that average.

---

## 15. Bottom line

The current project thesis can be sharpened from

```text
there seems to be an H4-like N^-13/8 correction
```

to a much more ambitious structural program:

```text
finite-size correction sectors may form representations of the exact Gaussian-cover
semigroup of the square CM torus; ordinary CFT fields are multiplicative eigenvectors,
LCFT fields are generalized eigenvectors with logarithmic norm cocycles, and matching
exchange supplies an independent block grading.
```

If that picture is right, the next breakthroughs should not come from ever finer free-exponent fits. They should come from **composition laws, cocycle ratios, modular-shape fingerprints, and tunable amplitude zeros**.

That is the level at which I would now try to break the problem.
