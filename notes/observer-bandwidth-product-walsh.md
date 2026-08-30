# Exact product-measure observer bandwidth

For independent Bernoulli sites with parameter `p`, write

`psi_S(x) = product_{i in S} (x_i-p)`.

These unnormalized centered products are mutually orthogonal and satisfy

`T_rho psi_S = rho^|S| psi_S`

exactly under retain-with-probability-`rho`, otherwise-resample-from-`p`
noise. Dividing by `[p(1-p)]^(|S|/2)` gives the normalized p-biased Walsh
basis. Therefore the covariance of a degree-`d` multilinear observer with any
square-integrable source is a polynomial in `rho` with nonconstant powers at
most `d`.

The committed rational oracle exhausts all 32 configurations for `n=5`, uses
`p=2/5`, and checks every basis function at three rational noise levels. That
is 3,072 exact pointwise eigenfunction checks. A synthetic degree-four
observer is then compared with an arbitrary full-cube source both by direct
transition-matrix summation and by its Walsh degree coefficients. All three
covariance residuals are exactly zero, and no coefficient above degree four
survives.

This is an observer/clock theorem only. It does not bound Gaussian-cover or
context-Hankel rank, identify continuum states, or turn the number of active
finite noise modes into Jordan evidence.
