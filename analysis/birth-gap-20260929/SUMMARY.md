# Executed joint-birth pilot — 2026-09-29

128,000 independent-filtration samples; 8 batches × 2,000 per lattice/size.
Local ARM64 acquisition including compilation, exact controls and benchmarks: 6.89 s, at most 4 concurrent jobs.

## Scale and direct jump

Errors after ± are complete-batch jackknife SE, not confidence limits.

| Lattice | L | Direct 0→2 / 16,000 | E[D]/L^(5/4) | E[D]/mixture IQR |
|---|---:|---:|---:|---:|
| square | 16 | 97 | 0.4230 ± 0.0032 | 0.8460 ± 0.0635 |
| square | 32 | 24 | 0.4206 ± 0.0029 | 0.7808 ± 0.0054 |
| square | 64 | 8 | 0.4273 ± 0.0026 | 0.7437 ± 0.0068 |
| square | 128 | 1 | 0.4278 ± 0.0030 | 0.7195 ± 0.0079 |
| triangular | 16 | 38 | 0.4032 ± 0.0018 | 0.8064 ± 0.0036 |
| triangular | 32 | 13 | 0.4031 ± 0.0027 | 0.7482 ± 0.0050 |
| triangular | 64 | 4 | 0.4075 ± 0.0030 | 0.7233 ± 0.0091 |
| triangular | 128 | 1 | 0.4056 ± 0.0031 | 0.7013 ± 0.0076 |

## Continuous-label near diagonal: P[L^(3/4)(T2−T1) ≤ ε]

W is the pooled J1/J2 interquartile width; zero atom included. Subtract the preceding direct-atom frequency for strictly positive mass.

| Lattice | L | ε=.025 | ε=.05 | ε=.1 | ε=.2 |
|---|---:|---:|---:|---:|---:|
| square | 16 | 0.05352 ± 0.00094 | 0.09853 ± 0.00139 | 0.18323 ± 0.00218 | 0.33486 ± 0.00290 |
| square | 32 | 0.04713 ± 0.00141 | 0.09132 ± 0.00193 | 0.17605 ± 0.00300 | 0.33208 ± 0.00319 |
| square | 64 | 0.04568 ± 0.00177 | 0.08830 ± 0.00320 | 0.17245 ± 0.00407 | 0.32576 ± 0.00314 |
| square | 128 | 0.04557 ± 0.00147 | 0.09160 ± 0.00148 | 0.17613 ± 0.00277 | 0.33549 ± 0.00193 |
| triangular | 16 | 0.05033 ± 0.00114 | 0.09686 ± 0.00185 | 0.18530 ± 0.00228 | 0.34564 ± 0.00129 |
| triangular | 32 | 0.04878 ± 0.00092 | 0.09536 ± 0.00172 | 0.18421 ± 0.00289 | 0.34440 ± 0.00411 |
| triangular | 64 | 0.04785 ± 0.00173 | 0.09253 ± 0.00203 | 0.17919 ± 0.00144 | 0.34160 ± 0.00341 |
| triangular | 128 | 0.04813 ± 0.00185 | 0.09363 ± 0.00218 | 0.17973 ± 0.00291 | 0.33807 ± 0.00390 |

## Exponent-free near diagonal: P[(N+1)(T2−T1)/W ≤ ε]

W is the pooled J1/J2 interquartile width; zero atom included. Subtract the preceding direct-atom frequency for strictly positive mass.

| Lattice | L | ε=.025 | ε=.05 | ε=.1 | ε=.2 |
|---|---:|---:|---:|---:|---:|
| square | 16 | 0.03009 ± 0.00220 | 0.05334 ± 0.00409 | 0.09818 ± 0.00764 | 0.18259 ± 0.01393 |
| square | 32 | 0.02623 ± 0.00097 | 0.05056 ± 0.00147 | 0.09796 ± 0.00201 | 0.18862 ± 0.00308 |
| square | 64 | 0.02707 ± 0.00111 | 0.05206 ± 0.00235 | 0.10100 ± 0.00396 | 0.19642 ± 0.00449 |
| square | 128 | 0.02728 ± 0.00104 | 0.05426 ± 0.00133 | 0.10843 ± 0.00175 | 0.20801 ± 0.00329 |
| triangular | 16 | 0.02639 ± 0.00066 | 0.05015 ± 0.00114 | 0.09650 ± 0.00184 | 0.18463 ± 0.00228 |
| triangular | 32 | 0.02669 ± 0.00053 | 0.05240 ± 0.00099 | 0.10232 ± 0.00181 | 0.19737 ± 0.00304 |
| triangular | 64 | 0.02721 ± 0.00140 | 0.05365 ± 0.00206 | 0.10369 ± 0.00269 | 0.20054 ± 0.00319 |
| triangular | 128 | 0.02832 ± 0.00154 | 0.05534 ± 0.00189 | 0.10771 ± 0.00262 | 0.20563 ± 0.00336 |

## Discrete near diagonal: P[D/L^(5/4) ≤ ε]

W is the pooled J1/J2 interquartile width; zero atom included. Subtract the preceding direct-atom frequency for strictly positive mass.

| Lattice | L | ε=.025 | ε=.05 | ε=.1 | ε=.2 |
|---|---:|---:|---:|---:|---:|
| square | 16 | 0.00606 ± 0.00058 | 0.06763 ± 0.00163 | 0.17763 ± 0.00303 | 0.32125 ± 0.00383 |
| square | 32 | 0.02569 ± 0.00118 | 0.07325 ± 0.00222 | 0.16375 ± 0.00368 | 0.33031 ± 0.00306 |
| square | 64 | 0.04106 ± 0.00159 | 0.08744 ± 0.00349 | 0.17250 ± 0.00468 | 0.32512 ± 0.00328 |
| square | 128 | 0.04200 ± 0.00160 | 0.08981 ± 0.00144 | 0.17588 ± 0.00338 | 0.33600 ± 0.00188 |
| triangular | 16 | 0.00237 ± 0.00039 | 0.06275 ± 0.00185 | 0.17694 ± 0.00337 | 0.33356 ± 0.00191 |
| triangular | 32 | 0.02569 ± 0.00081 | 0.07575 ± 0.00174 | 0.16975 ± 0.00276 | 0.34431 ± 0.00464 |
| triangular | 64 | 0.04306 ± 0.00193 | 0.09131 ± 0.00236 | 0.17969 ± 0.00169 | 0.34125 ± 0.00360 |
| triangular | 128 | 0.04519 ± 0.00187 | 0.09119 ± 0.00225 | 0.17913 ± 0.00305 | 0.33881 ± 0.00403 |

## Discrete exponent-free near diagonal: P[D/W ≤ ε]

W is the pooled J1/J2 interquartile width; zero atom included. Subtract the preceding direct-atom frequency for strictly positive mass.

| Lattice | L | ε=.025 | ε=.05 | ε=.1 | ε=.2 |
|---|---:|---:|---:|---:|---:|
| square | 16 | 0.00606 ± 0.00058 | 0.00606 ± 0.00058 | 0.06763 ± 0.00163 | 0.17763 ± 0.00303 |
| square | 32 | 0.02569 ± 0.00118 | 0.05025 ± 0.00203 | 0.09575 ± 0.00169 | 0.18550 ± 0.00393 |
| square | 64 | 0.02100 ± 0.00069 | 0.05100 ± 0.00207 | 0.09725 ± 0.00364 | 0.19112 ± 0.00975 |
| square | 128 | 0.02538 ± 0.00143 | 0.05025 ± 0.00173 | 0.10606 ± 0.00152 | 0.20738 ± 0.00282 |
| triangular | 16 | 0.00237 ± 0.00039 | 0.00237 ± 0.00039 | 0.06275 ± 0.00185 | 0.17694 ± 0.00337 |
| triangular | 32 | 0.02569 ± 0.00081 | 0.05250 ± 0.00130 | 0.10150 ± 0.00247 | 0.19369 ± 0.00377 |
| triangular | 64 | 0.02138 ± 0.00133 | 0.05306 ± 0.00202 | 0.10138 ± 0.00249 | 0.19631 ± 0.00157 |
| triangular | 128 | 0.02706 ± 0.00165 | 0.05381 ± 0.00193 | 0.10394 ± 0.00217 | 0.20287 ± 0.00308 |

## Soft near-diagonal Z(δ), G=L^(3/4)(T2−T1)

Z(δ)=E[(1−G/δ)+]; same samples and covariance, not additional independent evidence.

| Lattice | L | δ=.025 | δ=.05 | δ=.1 | δ=.2 |
|---|---:|---:|---:|---:|---:|
| square | 16 | 0.03005 ± 0.00070 | 0.05312 ± 0.00091 | 0.09728 ± 0.00132 | 0.17906 ± 0.00189 |
| square | 32 | 0.02443 ± 0.00085 | 0.04688 ± 0.00126 | 0.09047 ± 0.00183 | 0.17308 ± 0.00244 |
| square | 64 | 0.02356 ± 0.00084 | 0.04528 ± 0.00169 | 0.08800 ± 0.00265 | 0.16934 ± 0.00300 |
| square | 128 | 0.02288 ± 0.00091 | 0.04576 ± 0.00110 | 0.09013 ± 0.00152 | 0.17388 ± 0.00168 |
| triangular | 16 | 0.02644 ± 0.00065 | 0.05008 ± 0.00107 | 0.09582 ± 0.00161 | 0.18157 ± 0.00169 |
| triangular | 32 | 0.02482 ± 0.00047 | 0.04852 ± 0.00085 | 0.09434 ± 0.00151 | 0.18052 ± 0.00249 |
| triangular | 64 | 0.02409 ± 0.00101 | 0.04722 ± 0.00147 | 0.09170 ± 0.00146 | 0.17676 ± 0.00159 |
| triangular | 128 | 0.02444 ± 0.00128 | 0.04770 ± 0.00162 | 0.09254 ± 0.00209 | 0.17648 ± 0.00231 |

## Soft exponent-free Z(δ), G=(N+1)(T2−T1)/W

Z(δ)=E[(1−G/δ)+]; same samples and covariance, not additional independent evidence.

| Lattice | L | δ=.025 | δ=.05 | δ=.1 | δ=.2 |
|---|---:|---:|---:|---:|---:|
| square | 16 | 0.01815 ± 0.00123 | 0.02996 ± 0.00217 | 0.05294 ± 0.00402 | 0.09694 ± 0.00744 |
| square | 32 | 0.01389 ± 0.00055 | 0.02616 ± 0.00089 | 0.05027 ± 0.00131 | 0.09699 ± 0.00191 |
| square | 64 | 0.01402 ± 0.00052 | 0.02684 ± 0.00111 | 0.05168 ± 0.00219 | 0.10049 ± 0.00322 |
| square | 128 | 0.01365 ± 0.00054 | 0.02718 ± 0.00080 | 0.05440 ± 0.00102 | 0.10634 ± 0.00186 |
| triangular | 16 | 0.01439 ± 0.00043 | 0.02634 ± 0.00065 | 0.04990 ± 0.00107 | 0.09547 ± 0.00160 |
| triangular | 32 | 0.01372 ± 0.00031 | 0.02665 ± 0.00049 | 0.05209 ± 0.00090 | 0.10119 ± 0.00160 |
| triangular | 64 | 0.01360 ± 0.00079 | 0.02709 ± 0.00127 | 0.05294 ± 0.00181 | 0.10274 ± 0.00223 |
| triangular | 128 | 0.01435 ± 0.00103 | 0.02813 ± 0.00137 | 0.05487 ± 0.00174 | 0.10610 ± 0.00236 |

## Exact control and uncertainty

L=3 was enumerated once per lattice, all 9! permutations. Square: P(D=0)=3/35, E[D]=3/2, E[D²]=43/14. Triangular: 2/35, 3/2, 81/28. Every exact identity passed before the pilot.

JSON preserves full delete-one covariance, all leave-one estimates, atom Wilson intervals and an exact one-sided 95% upper bound if an atom count is zero. At n=16,000, zero observations would imply only P(atom)<0.0001872 at one-sided 95%; not zero probability. In the actual L=128 runs each lattice had one direct atom.

The Beta CDF uses P[Beta(d,N+1−d)≤x]=P[Bin(N,x)≥d]. Backward binomial-tail sums avoid subtraction near zero; d=0 is retained as an atom. This integrates over labels without adding Monte Carlo noise or evidence.

## What this does and does not decide

Read the accompanying research note for the interpretation. This pilot distinguishes shrinking single-step double births from finite near-diagonal occupancy. It does not establish absence of a diagonal atom in any scaling limit.
