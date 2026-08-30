# Preregistration: N680 same-lineage heldout

Frozen before generating or inspecting N680 data. N680 `(22+14i,26+2i)` is the exact next `1+i` child of N340; H4 flips negative again and the projective scalar remains a zero-amplitude control.

Fixed H4-amplitude forecasts:

- `two_mode_recurrence`: `-0.00184129704223` (target SE `0.000955`)
- `single_frozen_lambda0`: `-0.00102715380307` (target SE `8.57e-05`)
- `single_free_lambda`: `-0.00303781939565` (target SE `0.00102`)
- `scale_neutral`: `-0.00891196363712` (target SE `0.00076`)

N340 variance implies 3-sigma two-mode/free-single separation at 117659531 samples/shape. The frozen 120M design projects SE `0.0003949` and separations `3.030` (two/free), `2.061` (two/fixed), and `17.903` (two/neutral).

Power refers to fixed point forecasts. Source-fit uncertainty remains in a separate predictive score; no model, exponent, or basis may be refit after reveal.
