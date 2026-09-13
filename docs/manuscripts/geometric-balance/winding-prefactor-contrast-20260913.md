# The winding-prefactor contrast at widths 4, 8 and 12

2026-09-13. Computed for issue #741, returning to the #739 manuscript. This is a
numerical counterpart to #740, not a new mechanism programme and not a width scan.

## What was computed

The once-per-COMPLETE-component cylinder density

    nu_w^G(p) = expected number of horizontally winding components retired per vertical row

for two fixed inputs, both rigorously subcritical by elementary path bounds:

| graph | p | why subcritical |
|---|---|---|
| NN square site | 1/4 | `3p < 1` |
| matching NN+NNN site | 1/8 | `7p < 1` |

at `w = 4, 8, 12`. The contrast is

    R_4 = nu_4 * nu_12 / nu_8^2,      beta_eff(4) = log(R_4) / log(4/3).

If `nu_w = A w^-beta exp(-kappa w)(1+o(1))` then the exponential mass and the
amplitude cancel in `R_4` and `R_4 -> (4/3)^beta`. This is therefore a direct
diagnostic of the `w^{-1/2}` sewing hypothesis of #740 — one number, not a fit of
three unknown parameters.

## Result

| graph | p | log nu_4 | log nu_8 | log nu_12 | beta_eff | deviation from 1/2 |
|---|---|---|---|---|---|---|
| NN | 1/4 | −5.419513505388 | −10.023930900251 | −14.400335955877 | **0.792584457** ± 1.9e−9 | +0.292584 |
| matching NN+NNN | 1/8 | −5.742669143987 | −10.180101524507 | −14.471261214034 | **0.508452577** ± 5.1e−10 | +0.008453 |

`w = 4` and `w = 8` are exact rationals with a zero certificate bound; `w = 12`
carries the certificate below, all below 6e−10 on `log nu`.

### The three-width agreement is not an asymptotic proof, and here it is not even self-consistent

Two checks make that concrete rather than rhetorical.

1. **A second, equally admissible window disagrees wildly.** The same construction on
   `(2,4,8)` returns `beta_eff = −7.137` (NN) and `−7.070` (matching). The assumed
   form does not describe those widths at all, which says the `(4,8,12)` values are
   finite-window effective exponents rather than a settled amplitude.
2. **The effective kappa is still moving.** Adjacent differences give

       NN        (4,8) 1.151104   (8,12) 1.094101      matching  (4,8) 1.109358   (8,12) 1.072790

   both still decreasing in magnitude at `w = 12`, so the amplitude has not settled either.

An identical construction cannot have two different true exponents, so at least one
of the two `beta_eff` values above is not asymptotic. The honest reading is:

- the `1/2` sewing hypothesis is **numerically close** on the matching graph at `p = 1/8`;
- it is **not** close on NN at `p = 1/4`, where the same formula returns `0.79`;
- and the `(2,4,8)` window shows the whole two-parameter description is not yet valid,
  so neither number should be quoted as `beta`.

Reported as required even though inconsistent with `1/2`. This neither confirms nor
refutes the sewing hypothesis; deciding that needs the `A, beta` derivation of #740,
not a fourth width.

## Engine, validation and cost

The supplied `scripts/cylinder_winding_intensity.py` caps the builder at width ≤ 10 and
solves densely in `Fraction`. `scripts/cylinder_winding_intensity_fast.cpp` is an
allocation-free C++ port of the same `advance()` / `empty_state()` / `reward_lump()`
semantics, with the transition table cached after BFS so that the reward-lumping
refinement does not recompute 6e8 successors per iteration.

**Validation is exact, not indicative.** All 18 published controls of
`results/geometric-consistency/cylinder-winding-intensity.json` — both graphs, widths
2/3/4, `p = 1/4, 1/2, 3/4` — reproduce as **equal rationals**, and the frontier-state
and reward-lump counts match (`6/3`, `14/4`, `38/7`). The capacity-probe counts
`102/282/786/2214` at `w = 5..8` and `90` reward lumps at `w = 8` also reproduce.

One caveat worth recording, found by this validation: the C++ port first disagreed on
the matching graph alone. The cause was integer division — C++ `/` truncates toward
zero while Python `//` floors, and the only affected call is the `dx = −1` diagonal step
at `i = 0`, where `(i+dx)//w = −1` but `(i+dx)/w = 0`. NN was unaffected, which is
exactly why a matching-side control matters. Fixed by an explicit floor division.

Cost, measured rather than extrapolated, single process on a 16 vCPU aarch64 container:

| | w = 4 | w = 8 | w = 12 |
|---|---:|---:|---:|
| frontier states | 38 | 2 214 | 147 578 |
| reward lumps | 7 | 90 | 2 105 |
| transitions | 152 | 566 784 | 604 479 488 |
| wall time | <1 s | 0.5 s | 504 s (NN) / 561 s (matching) |
| peak RSS | — | — | 2.4 GB (2.25 GB cached table) |

`w = 12` broke down as BFS 234 s + table cache 238 s + refinement ≈ 30 s. The same
counts are reachable on a laptop at `w = 8` (0.5 s) but not at `w = 12`.

## Error control

For a normalised nonnegative `pi_hat` and the uniform empty-row reset
`delta = (1-p)^w` (the empty row always maps to the empty state, so the chain is
uniformly ergodic with gap at least `delta`),

    |pi_hat . g - nu| <= ||g||_inf * ||pi_hat K - pi_hat||_1 / delta.

`w = 4, 8` are solved in exact rational arithmetic. `w = 12` is solved in float64
(2 105 states), then corrected once using the **exact rational** stationary residual and
certified exactly; that step is what brings the bound from 1.6e−14 down to 3e−16 (NN)
and 7.5e−17 (matching), i.e. `log nu` to 5.5e−10 and 1.5e−10.

## What this does not establish

No exponent is determined. No asymptote is claimed. The `w^-beta` form is not verified;
the `(2,4,8)` window shows it fails at small width. `kappa_G` is not measured as a
limit. No new `p_c`, no Monte Carlo, no GPU, no continuum identification. The residual
certificate bounds arithmetic only; it says nothing about the model assumptions or the
imported inputs of #739.

## Reproducing

```sh
g++ -O3 -std=c++17 -o winding_build scripts/cylinder_winding_intensity_fast.cpp
for w in 4 8 12; do ./winding_build $w 0 out_nn_$w.json; done    # NN, 0 = site NN
for w in 4 8 12; do ./winding_build $w 1 out_m_$w.json;  done    # matching
python3 scripts/winding_prefactor_contrast.py \
    '[["out_nn_4.json","1/4",true],["out_nn_8.json","1/4",true],["out_nn_12.json","1/4",false],'\
'"["out_m_4.json","1/8",true],["out_m_8.json","1/8",true],["out_m_12.json","1/8",false]]' \
    results/geometric-consistency/winding-prefactor-contrast.json
```

The `w = 12` builds need ~2.4 GB and ~9 minutes each; the `w ≤ 8` builds are seconds.
