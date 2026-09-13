# E_top: modulus interaction versus direction transport

## Geometry changes the response direction, not only its amplitude

The completed N50 archive rejects one common amplitude for the four coordinates A_top,E_top,C,W. This is a new retrospective model comparison using zero new samples; it is not an E_top-only detection or field identification.

Free common amplitude: lambda=3.32170792; chi2=19.8718/3, p=0.000180452.
The lambda profile interval belongs to the rejected model and is not used in the predictions below.

A and E are fixed-p rank probabilities; C and W are integrated clock/lifetime coordinates. Every input is a first-minus-second direction contrast divided by its exact cos(4 theta) contrast, at N50 and p_ref=.59274605079. Two random streams are independent, while each four-vector retains its full covariance.

| Declared common-vector shape | lambda | chi2/4 | p |
|---|---:|---:|---:|
| no_modulus_response | 1 | 236.756 | 4.63371e-50 |
| height | 2 | 89.4961 | 1.68479e-18 |
| height_squared | 4 | 37.5505 | 1.38716e-07 |
| pure_area_E4 | 2.75 | 32.7538 | 1.34142e-06 |

## Is a second transport amplitude needed?

Odd A/C lambda=3.488293; even E/W lambda=7.522162.
Two-amplitude chi2=8.7014/2, p=0.0128978; improvement=11.1704/1, p=0.000831142.

The parity split is a narrow survivor at alpha=.01 and fails at .05. It is a useful direction, not a selected mechanism. Earlier A+C-plane compatibility is not contradicted: a four-vector can change direction while still satisfying E=beta*A+gamma*C. That earlier test used different geometries/centers and is not pooled here.

## A sheared geometry separates the conditional shape hypotheses

A vector-affine family uses v(tau)=(1-t)x_i+t*x_2i with independent offset and slope vectors, so the source vectors need not share a ray. The four choices of t are log(height), height squared, height-only E4, and the actual complex-modulus E4 covector. Each exactly interpolates the two source means; none is selected by those two points.

For N100, assume both affine coefficients for coordinate j receive the same unknown gain c_j(N100). This is an explicit new shape/area-separability hypothesis. In the SAME-N100 ratio to a new 2i bridge, that gain cancels:

`R_j(tau)/R_j(2i) = t+(1-t)*x_j(50,i)/x_j(50,2i)`.

No N^-13/8 or other area exponent is imposed. The old-N50 affine vectors in JSON are prototypes, not absolute N100 amplitude predictions. Source-only 95% Fieller sets below retain weak-denominator nonidentification; they are coordinatewise, not simultaneous.

| Shape family | target | coordinate | same-N100 ratio | source 95% Fieller set |
|---|---|---|---:|---|
| affine_log_height | tau_4i | A_top | 1.679091 | [1.41495, 1.78523] |
| affine_log_height | tau_4i | E_top | 1.567161 | [-inf, 1.83381] U [3.07301, +inf] |
| affine_log_height | tau_4i | C | 1.705493 | [1.64814, 1.74755] |
| affine_log_height | tau_4i | W | 1.807015 | [1.63049, 1.87368] |
| affine_log_height | tau_half_plus_i | A_top | 0.3209095 | [0.214767, 0.585047] |
| affine_log_height | tau_half_plus_i | E_top | 0.4328392 | [-inf, -1.07301] U [0.166194, +inf] |
| affine_log_height | tau_half_plus_i | C | 0.2945068 | [0.252449, 0.351858] |
| affine_log_height | tau_half_plus_i | W | 0.1929845 | [0.126319, 0.369507] |
| affine_height_squared | tau_4i | A_top | 3.716362 | [2.65981, 4.14093] |
| affine_height_squared | tau_4i | E_top | 3.268643 | [-inf, 4.33522] U [9.29205, +inf] |
| affine_height_squared | tau_4i | C | 3.821973 | [3.59257, 3.99021] |
| affine_height_squared | tau_4i | W | 4.228062 | [3.52197, 4.49472] |
| affine_height_squared | tau_half_plus_i | A_top | 0.3209095 | [0.214767, 0.585047] |
| affine_height_squared | tau_half_plus_i | E_top | 0.4328392 | [-inf, -1.07301] U [0.166194, +inf] |
| affine_height_squared | tau_half_plus_i | C | 0.2945068 | [0.252449, 0.351858] |
| affine_height_squared | tau_half_plus_i | W | 0.1929845 | [0.126319, 0.369507] |
| affine_height_E4 | tau_4i | A_top | 4.197857 | [2.95402, 4.69769] |
| affine_height_E4 | tau_4i | E_top | 3.670777 | [-inf, 4.92642] U [10.7619, +inf] |
| affine_height_E4 | tau_4i | C | 4.322188 | [4.05212, 4.52024] |
| affine_height_E4 | tau_4i | W | 4.800259 | [3.96901, 5.11419] |
| affine_height_E4 | tau_half_plus_i | A_top | 0.3209095 | [0.214767, 0.585047] |
| affine_height_E4 | tau_half_plus_i | E_top | 0.4328392 | [-inf, -1.07301] U [0.166194, +inf] |
| affine_height_E4 | tau_half_plus_i | C | 0.2945068 | [0.252449, 0.351858] |
| affine_height_E4 | tau_half_plus_i | W | 0.1929845 | [0.126319, 0.369507] |
| affine_E4 | tau_4i | A_top | 4.197857 | [2.95402, 4.69769] |
| affine_E4 | tau_4i | E_top | 3.670777 | [-inf, 4.92642] U [10.7619, +inf] |
| affine_E4 | tau_4i | C | 4.322188 | [4.05212, 4.52024] |
| affine_E4 | tau_4i | W | 4.800259 | [3.96901, 5.11419] |
| affine_E4 | tau_half_plus_i | A_top | 0.08194691 | [-0.0615458, 0.439031] |
| affine_E4 | tau_half_plus_i | E_top | 0.2332631 | [-inf, -1.80248] U [-0.127211, +inf] |
| affine_E4 | tau_half_plus_i | C | 0.04625342 | [-0.0106045, 0.123786] |
| affine_E4 | tau_half_plus_i | W | -0.09099306 | [-0.181117, 0.147646] |

The cleanest source-resolved contrast is C: actual-modulus affine-E4 predicts shear/2i=0.0462534, whereas every listed height-only family predicts 0.294507. Their source uncertainty sets are distinct; future target noise and model discrepancy remain to be measured. E_top's denominator is weak, so its unbounded Fieller set cannot be presented as a precise forecast. A possible W sign reversal is a hypothesis, not a resolved sign prediction.

N50 with Smith(1,50)/(5,10) admits no third modulus. N100 with Smith(1,100)/(5,20) admits 2i,4i,1/2+i. The exact proof/matrices are at commit b9e4ea1. The bridge-and-shear pair is a smaller source-informed test; all three N100 shapes permit a new-area offset/amplitude-free test that drops the cross-area separability assumption.

## Dependence and scope

The common four-vector amplitude is rejected. Its lambda interval describes that rejected family and is NOT used for forecasts. Two source shapes saturate each vector-affine family, so none is source-selected. Pure E4 common-ray rejection does not exclude all thermal Q4 couplings. Common/parity scores use approximate Gaussian covariance inference.

The likelihood profiles the latent source vector with both 4x4 covariance blocks. It does not treat the source point as exact, add a residual log-determinant, combine prior p-values, or rerun the factorial score. The newly selected model comparisons and shape forecasts are retrospective, exploratory outputs, not preregistered source evidence.

Source: `087b07cb69f2481cbcd55fe2194150d7620835e5:results/p267-etop-tau-topology-factorial/score.json`, SHA256 `379eedf05559f259a237b10bc0a5b0d3e26540cd78f222218d4d31a56ae4467a`.

Reproduce with `/Users/lc/python-envs/research-py311/bin/python scripts/etop_modulus_survivors.py`. The pinned source Git object must be present (fetch the PR267 source branch if needed).
