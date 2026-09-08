# Two-lobe Gaussian moment closure

Post-reveal deterministic mechanism model chosen after observing two peaks and standardized moments; no independent new evidence, no Gaussian assumption imposed on signed data. Moments 3/4 identify two lobe parameters and moments 5/6 are the overidentifying readout. Source fitting repeated inside every aligned delete-one replicate.

| Quantity | N100 | N400 | change | change SE |
|---|---:|---:|---:|---:|
| right_lobe_weight | 0.380789591 | 0.435950881 | 0.05516129 | 0.0207348 |
| between_lobes_variance_fraction | 0.883463277 | 0.858490446 | -0.024972831 | 0.0142307 |
| within_lobe_sigma | 0.341374755 | 0.376177556 | 0.034802802 | 0.019079 |
| left_lobe_center | -0.737085149 | -0.814569323 | -0.077484175 | 0.0349263 |
| right_lobe_center | 1.19859053 | 1.05391944 | -0.14467109 | 0.0461899 |
| predicted_mu5 | 1.28240357 | 0.65535003 | -0.62705354 | 0.229977 |
| predicted_mu6 | 3.48528832 | 3.24710265 | -0.23818567 | 0.300306 |
| mu5_residual | -0.153414942 | -0.170387122 | -0.01697218 | 0.0195704 |
| mu6_residual | -0.198722616 | -0.0100164819 | 0.18870613 | 0.0483648 |

N100 residual moments 5/6: {'chi2': 848.3869780844277, 'df': 2, 'nominal_p': 5.958109076448906e-185}

N400 residual moments 5/6: {'chi2': 177.59403345921004, 'df': 2, 'nominal_p': 2.7286361269471746e-39}

Cross-scale residual change: {'chi2': 46.742700107670615, 'df': 2, 'nominal_p': 7.078669519694668e-11}

## Any common symmetric lobe kernel

Structural extension after the Gaussian closure result: moments 3/5 identify the two-center component for any common symmetric lobe kernel. Its reconstructed sixth even moment and 2x2 moment Hankel determinant must be nonnegative. No extra fitted kernel family.

| Derived necessary quantity | N100 | SE | N400 | SE |
|---|---:|---:|---:|---:|
| right_lobe_weight | 0.388833784 | 0.00727445 | 0.445767727 | 0.0185376 |
| between_lobes_variance_fraction | 0.928003062 | 0.00425629 | 0.96069976 | 0.0272568 |
| kernel_variance | 0.071996938 | 0.00425629 | 0.0393002397 | 0.0272568 |
| kernel_fourth_moment | 0.185946422 | 0.0066983 | 0.381742349 | 0.0924043 |
| kernel_sixth_moment | -1.75780018 | 0.0741887 | -3.84934607 | 1.20234 |
| kernel_moment_Hankel_determinant | -0.161132302 | 0.0081322 | -0.297007444 | 0.0507798 |
