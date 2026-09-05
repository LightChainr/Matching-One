# New N100 thermal deformation: translation, dilation, or another shape?

The lowest odd signed moments determine a tangent-transport velocity. The remaining odd moments and the entire even profile are then genuine additional model constraints, with all same-stream source uncertainty retained. This is post-reveal mechanism analysis, not fresh evidence or a test-suite rerun.

## Odd profile: low moments fix the transport, higher moments constrain it

| candidate | anchored moments | theta | remaining chi-square / df | nominal p |
|---|---|---|---:|---:|
| translation_tangent | [1] | [0.07188292008613974] | 481.5019 / 5 | 7.84531e-102 |
| affine_velocity_tangent | [1, 2] | [0.05033264006899235, -0.051059950392760796] | 217.8451 / 4 | 5.45283e-46 |
| quadratic_velocity_tangent | [1, 2, 3] | [0.12690858698444896, -0.0725144927042187, -0.05291981071087471] | 148.2432 / 3 | 6.30551e-32 |

## Even profile: transfer the same normalized velocity, not a new fit

Each parity has a separate area amplitude r=U_0/D_0. Odd anchors determine theta/r_A; the even moments use r_E times this same velocity.

| candidate | even remaining chi-square / df | nominal p |
|---|---:|---:|
| translation_tangent | 795.9493 / 6 | 1.15525e-168 |
| affine_velocity_tangent | 295.5528 / 6 | 7.33791e-61 |
| quadratic_velocity_tangent | 156.3153 / 6 | 3.56943e-31 |

All uncertainty is estimated by deleting each of the 200 aligned batches and refitting its area amplitudes and low-moment velocity. Gaussian-reference scores are exploratory; the candidates are not independent evidence blocks.

## What the model means

Write D=Y(4i)-Y(2i), U=Y(1/2+i)-Y(2i). After fixing the area ratio r, the candidate is U-rD=(1/s) d_p[v(z)D], with z=s(p-p_ref). Integration by parts gives R_j=-j sum_l theta_l D_(j-1+l). Constant v is a translation tangent; linear v adds dilation with the area-preserving amplitude correction; quadratic v introduces a non-affine tangent.

This is a signed-profile transport model, not an assumption that A_top is a probability density. Ordinary observable reparameterization omits the Jacobian. The present scores do not settle finite nonlinear monotone transport, and cannot identify a continuous field.
