# Erratum — the (2,4,8) window and the κ-drift reading in `winding-prefactor-contrast.json`

**Status: erratum to `results/geometric-consistency/winding-prefactor-contrast.json` (commit `3745b13`)
and to the interpretation printed in the `3745b13` commit message, the #741 delivery comment, and the
matching #739 comment. The committed JSON and note are left untouched so the audit trail shows what
was actually computed; this file is the correction of record. Source of the correction: #741 comment
5652002027. Both corrections were re-verified here against the committed log values before being
accepted.**

## E1. The (2,4,8) window used an equal-spacing formula on unequal spacing

The committed fields `beta_eff_2_4_8` (−7.137336 NN, −7.069824 matching) and `R_window_2_4_8` were
computed as `log nu_2 − 2 log nu_4 + log nu_8`, the second difference that cancels the mass term
**only for equally spaced widths**. The window (2,4,8) has spacing 2 then 4, so the combination
retains −2κ: it is not a mass-cancelled contrast, and the values are not beta estimates of any kind.
The conclusion drawn from them — "the assumed form does not describe those widths at all" — is
therefore **retracted**.

The correct contrast for distinct widths x < y < z uses c = (z−y, x−z, y−x), which satisfies both
Σc_i = 0 and Σc_i·w_i = 0, and divides by −Σc_i·log w_i:

```
beta_eff(x,y,z) = [ (z−y)·log nu_x + (x−z)·log nu_y + (y−x)·log nu_z ] / [ −Σ c_i log w_i ]
```

For (2,4,8) this is `(2·log nu_2 − 3·log nu_4 + log nu_8)/log 2`, giving

| graph | committed (wrong) | corrected |
|---|---|---|
| NN, p = 1/4 | −7.137336 | **0.718245787** |
| matching, p = 1/8 | −7.069824 | **0.533377434** |

Re-verified here at 40 digits against the committed `log_nu_by_width` values; agreement with the
corrector's numbers is exact to 1e−8.

The substantive reading changes accordingly: the two windows now read (2,4,8) → 0.718/0.533 and
(4,8,12) → 0.793/0.508. That is **finite-window drift between two valid mass-cancelled contrasts**,
not evidence that the two-parameter form breaks down, and not a second exponential. The (4,8,12)
values, the six densities and the engine validation in the committed JSON are unaffected (that window
IS equally spaced and the formula used there is correct).

## E2. The κ-drift reading was vacuous, not evidence about the amplitude

The `3745b13` comment said "the effective kappa is still moving at w = 12 … Neither graph has settled
its amplitude." But for an **exact** `A·w^−β·e^(−κw)` law, adjacent κ estimates must drift by

```
κ_eff(u,v) − κ_eff(v,z) = β·log(4/3)/4        (for u,v,z equally spaced with step u/2… here 4)
```

so the drift is forced by the fitted β itself and carries no information beyond `beta_eff`. Checked
against the committed values:

| graph | observed drift | β-forced prediction | difference |
|---|---|---|---|
| NN | 0.057003 | 0.057003 | −3.9e−16 |
| matching | 0.036568 | 0.036568 | −1.4e−16 |

The observed drifts equal the β-forced amounts to machine precision. The "amplitude has not settled"
reading is **retracted**; what remains is only the honest statement that three widths cannot
distinguish a true asymptote from a finite-window effective exponent, which the ticket already said.

## What survives unchanged

- The six exact/certified densities and their certificates.
- The (4,8,12) contrasts 0.792584457 (NN) / 0.508452577 (matching) and the observation that the two
  graphs disagree.
- The engine, the 18-control validation, the w=12 cost report, and the floordiv portability note.

Filed on the record by the same author who made the errors, after independent re-verification.
