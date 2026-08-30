# Preregistration: N170 exact angle-flip natural current

Frozen before generating or inspecting any N170 projective-birth data.

N170 `(11+7i,13+i)` is the exact `1+i` child of N85. The reflection-even H4 covector flips sign for both orientations, while the charged scalar `q_A^2=(u+H_F3)/2` remains fixed.

The H4-only model from commit `186d72a` freezes:

- `K_A(first)=+0.00393553360239`
- `K_A(second)=-0.00585020935024`
- `Delta_K_A(second-first)=-0.00978574295262`
- H4 amplitude `-0.00613732576673`
- charged/projective scalar `0`

The primary parameter-free discriminator is the sign-flipped negative pair contrast versus scalar/common-mode zero. Residual in the H4 amplitude coordinate is called curvature; residual in the orthogonal A-scalar coordinate is called projective/common mode.

Design: 8M samples/shape, 80 aligned batches, seed `202608337170`, on idle Huawei `DevEnvC_HZsCM6`. N145 covariance projects pair SE `0.001718` and expected H4-versus-zero `z=5.696`.

No H4/H8 vote, exponent fit, or post-reveal basis change is allowed.
