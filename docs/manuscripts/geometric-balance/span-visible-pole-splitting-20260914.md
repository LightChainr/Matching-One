# Coalescing visible poles in the complete-component span spectrum

Date: 2026-09-14

Status: bounded diagnostic from the already committed `span-spectrum-20260913.json`, plus a correction to one proposed proof route for #758/#760. No new transfer build or width production is used.

## 1. Why this check matters

For fixed circumference `w`, the tagged complete-component span sequence has a finite-state representation

```text
d_h = alpha R^(h-1) b.
```

Hence its ultimate `h->infinity` decay is controlled by the slowest **visible** pole of the source/readout pair. A tempting bridge between #758 and #760 is to assume that the relaxation to that pole occurs on a scale `o(w)`, so that the simultaneous macroscopic regime `h=A w` already sees the fixed-w Perron tail.

The existing span data give a strong negative diagnostic for that assumption.

## 2. Frozen two-pole extraction

Input:

```text
results/geometric-consistency/span-spectrum-20260913.json
```

No solver is rerun.

For each stored `d_h` sequence, fit on the late window beginning at `h=13` the second-order recurrence

```text
d_(h+2) = S d_(h+1) - P d_h.                              (2.1)
```

Each row is divided by `d_(h+1)` before the two-parameter least-squares solve, avoiding the many-decade dynamic range. The two recurrence roots are reported as effective visible poles

```text
rho_1 > rho_2 > 0.
```

A separate two-exponential fit

```text
d_h ~= c_1 rho_1^(h-1) + c_2 rho_2^(h-1)                 (2.2)
```

checks the residues. This is a projected/Prony diagnostic, not a certified full-matrix eigensolve.

For NN `p=1/4`, widths 2--7 give stable real roots and late-window relative residuals between about `1e-12` and `1e-7` (the largest widths inherit the stored float residuals).

## 3. NN p=1/4 result

Representative values are:

| w | rho1 | rho2 | Delta gamma=log(rho1/rho2) | xi_rel=1/Delta gamma | xi_rel/w | c2/c1 |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | .290359457 | .187500000 | .437340815 | 2.29 | 1.14 | -.9984 |
| 3 | .329288918 | .303004745 | .083187071 | 12.02 | 4.01 | -1.0049 |
| 4 | .344536629 | .334361548 | .029977523 | 33.36 | 8.34 | -1.0204 |
| 5 | .350913754 | .346766728 | .011888178 | 84.12 | 16.82 | -1.0153 |
| 6 | .353568560 | .351836957 | .004909537 | 203.69 | 33.95 | -1.0091 |
| 7 | .354597156 | .354034203 | .001588845 | 629.39 | 89.91 | -1.0038 |

Two features are much more important than the individual digits:

1. the two slow visible poles rapidly coalesce;
2. their fitted residues are nearly equal and opposite.

Thus the scalar span response is not approaching its one-pole regime on an `O(w)` height scale in these controls.

The same phenomenon becomes even stronger deeper in the subcritical phase / on the matching graph. Examples from the same stored file:

```text
NN p=1/8, w=6:        xi_rel ~= 5288,  c2/c1 ~= -1.00026;
matching p=1/16, w=5: xi_rel ~= 321,   c2/c1 ~= -1.00365;
matching p=1/16, w=6: xi_rel ~= 1209,  c2/c1 ~= -1.00138.
```

So this is not a one-point anomaly of NN `p=1/4`.

## 4. Direct macroscopic-hazard check

Let

```text
T_w(h)=sum_{j>=h} d_j
```

with the stored tail bin appended. Define the finite-height hazard

```text
gamma_eff,w(h) = -log[T_w(h+1)/T_w(h)].                   (4.1)
```

For NN `p=1/4`:

```text
w=4: gamma_eff(w)=0.7790, gamma_eff(2w)=0.9436,
     ultimate two-pole gamma1=1.0656;
w=5: 0.7867, 0.9388, gamma1=1.0472;
w=6: 0.8042, 0.9440, gamma1=1.0397;
w=7: 0.8242, 0.9519, gamma1=1.0368.
```

The `h=w,2w` tails are visibly not in the fixed-w ultimate Perron regime.

This does **not** falsify the #758 macroscopic LDP. It falsifies the shortcut that would derive that LDP by first taking `h->infinity` at fixed w and then assuming a relaxation length `o(w)`.

## 5. Noncommuting limits and a replacement asymptotic picture

Write the two slow pole masses as

```text
gamma_1 = gamma_bar - Delta/2,
gamma_2 = gamma_bar + Delta/2,
```

and suppose the corresponding residues are approximately `A` and `-A`, plus a smaller common component `B/2`. Then

```text
d_h
 ~= exp(-gamma_bar h)
    [ B cosh(Delta h/2) + 2 A sinh(Delta h/2) ].           (5.1)
```

If `Delta w -> 0`, then for `h=A0 w`

```text
d_h
 ~= exp(-gamma_bar h)
    [ B + A Delta h + o(Delta h) ].                       (5.2)
```

Thus a coalescing pair can have an exponentially long fixed-w relaxation time while the simultaneous `h=O(w)` exponential rate is still governed by the common centre mass `gamma_bar`. The near-cancellation changes subexponential/polynomial factors much more strongly than the large-deviation exponent.

This is a better interface for #758:

> prove that all macroscopically relevant visible slow modes have masses converging to the same planar cost `kappa`, while all truly fast modes remain separated; do **not** require the slow doublet to have mixed to one Perron vector by height `A w`.

## 6. Relation to #760

The periodic-locality argument on #760 concerns the location of the ultimate cylinder mass `gamma_w` relative to the plane mass `kappa`. The coalescing-pole diagnostic does not contradict an exponentially small `gamma_w-kappa` correction.

What it changes is measurement/proof strategy:

- finite-height hazards at `h=O(w)` should not be used as direct estimators of the ultimate fixed-w `gamma_w`;
- conversely, the fixed-w leading Perron pole should not be substituted into an `h=A w` morphology argument without controlling the neighboring visible pole and its residue.

The two limits

```text
h->infinity at fixed w,
and
w->infinity with h=A w
```

need not commute at the level of the visible decomposition even if both exponential masses converge to `kappa`.

## 7. A new metastable-doublet conjecture

For fixed subcritical p, the data suggest a two-dimensional slow visible subspace whose splitting is exponentially small in w. A deliberately weak conjecture is

```text
Delta gamma_w = exp[-Theta(w)],                            (7.1)
```

with opposite-sign source/readout residues and a common centre mass tending to the plane connection mass.

For NN `p=1/4`, the observed splitting is numerically comparable to `w nu_w`:

```text
Delta gamma_w/(w nu_w)
 = 1.83, 1.69, 1.77, 1.94, 1.67     for w=3,...,7.
```

At other p/graphs the proportionality constant is very different, so no universal constant is claimed. A more plausible structural version is

```text
Delta gamma_w = Theta_p(essential-component coverage/renewal intensity),
```

with `w nu_w` only a first proxy.

A possible mechanism is a metastable two-sector transfer block: rare complete-topology events weakly couple two otherwise degenerate longitudinal sectors, producing a symmetric/antisymmetric splitting and the observed opposite residues. The physical sectors should be identified from the actual tagged states before this interpretation is promoted.

## 8. What should be checked next

1. Extract the two leading **visible** poles directly from the tagged resolvent/operator, with left/right source overlaps, rather than only from `d_h`.
2. Identify which tagged frontier states carry the two slow eigenvectors.
3. Compare `Delta gamma_w` to `nu_w E_Palm K`, `nu_w E_Palm L`, and other renewal/coverage intensities once #774 supplies the relevant moments.
4. Repeat the projected recurrence on the stored matching `p=1/8` sequence and any future width-8 span spectrum without rebuilding the transfer.
5. For #758, formulate the upper-bound/LDP proof directly in the simultaneous `h=A w` limit; the coalescing pair explains why a fixed-w single-pole shortcut is unsafe.

## 9. Claim boundary

- Data statement: the two-pole recurrence and residue diagnostics on the committed span spectra.
- Strong negative diagnostic: `xi_rel=o(w)` is incompatible with the observed width trend and should not be used as a proof assumption without a new mechanism that reverses it.
- Conjecture: exponential slow-mode coalescence / metastable two-sector interpretation and its relation to renewal intensity.
- Not claimed: a certified full transfer eigendecomposition, an asymptotic exponent for the splitting, or a modification of the #758 rate function.
