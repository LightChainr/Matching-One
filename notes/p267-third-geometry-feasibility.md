# P267: a third modulus is impossible at N50, and the minimal replacement is N100

Source: `ce01e4d10abb03abf1b278192510937ee96db29d`, the completed N50
tau-by-topology-map factorial. This is a finite geometric result and a
conditional prediction basis, not new production or an operator identification.

## Exact finite-set result

Retain the map `O=(1/5)[[4,-3],[3,4]]` and the Smith pair
`(1,N) -> (5,N/5)`. Every target matrix is `P'=5B`, with primitive integer
`B` and `det B=N/25`. Conversely, every such B for which
`P=[[4,3],[-3,4]] B` has entry gcd one is an admissible pair.
Right-unimodular changes only change the period basis, so enumerating column
HNF B is exhaustive, not a matrix-box search.

| N | det B | complete modulus set, identifying reflection | number |
|---:|---:|---|---:|
| 25 | 1 | i | 1 |
| 50 | 2 | i, 2i | 2 |
| 75 | 3 | 3i, 1/2+3i/2 | 2 |
| 100 | 4 | 2i, 4i, 1/2+i | 3 |
| 125 | 5 | i, 5i, 1/2+5i/2 | 3 |

In particular, at N50 an index-two B has only three HNFs. Two give the
rectangle 2i and one the square i. **No third modulus exists with the
Smith-(5,10) endpoint, even if O is changed.** Cyclic-to-cyclic pairs would
change the original factorial question, not add its missing third cell.
N75 is the next admissible area but still has only two shapes. N100 is the
smallest three-shape design in the fixed-map Smith family. At N100 the
square endpoint is 10 times a rotation of the unit lattice: its Smith type
is (10,10), and its O-preimage has type (2,50), not the required pair.

## Runnable exact periods

Matrices are row-major; columns are period vectors. Every second matrix
below is exactly O times the first. Delta is first-minus-second cos(4 theta),
where theta is the first period's angle. All deltas are nonzero. All immediate
NN and matching neighbours are distinct; these cells have no local loop or
parallel-neighbour degeneration. The angle and Gram data are saved in
`results/p267-third-geometry-feasibility/certificate.json`.

| N | tau | cyclic P | noncyclic OP | delta chi4 | shortest period squared |
|---:|---|---|---|---|---:|
| 100 | 2i | `[[7,-2],[1,14]]` | `[[5,-10],[5,10]]` | 1152/625 | 50 |
| 100 | 4i | `[[4,12],[-3,16]]` | `[[5,0],[0,20]]` | -1152/625 | 25 |
| 100 | 1/2+i | `[[8,10],[-6,5]]` | `[[10,5],[0,10]]` | -1152/625 | 100 |
| 125 | i | `[[11,2],[-2,11]]` | `[[10,-5],[5,10]]` | 16128/15625 | 125 |
| 125 | 5i | `[[4,15],[-3,20]]` | `[[5,0],[0,25]]` | -1152/625 | 25 |
| 125 | 1/2+5i/2 | `[[7,1],[1,18]]` | `[[5,-10],[5,15]]` | 1152/625 | 50 |

The existing `src/threshold_rank_integer_period_mc.cpp` directly accepts each
row through `--first-matrix A B C D --second-matrix A B C D`. The certificate
stores these exact argument lists. No runner modification or execution is
part of this result.

## The allowed spin-four shape covector

Write complex periods omega1,omega2, tau=omega2/omega1 and y=Im(tau).
The basis-invariant covector is

`Phi4(P) = y^2 exp(-4i arg(omega1)) E4(tau)`.

This follows from `N^2 omega1^(-4) E4(tau)` and the weight-four modular
transformation. The normalized real map difference is

`Re[Phi4(P)-Phi4(OP)] / (chi4(P)-chi4(OP)) / E4(i)`.

For **these** cells Re(tau) is zero or one-half, so E4(tau) is real and this
reduces to `g(tau)=y^2 E4(tau)/E4(i)`. For a general sheared torus E4 is
complex: dropping its imaginary contribution would be incorrect. The script
refuses that unqualified real-only simplification.

| tau | g(tau) | t=(g-1)/(7/4) |
|---|---:|---:|
| i | 1 | 0 |
| 2i | 2.75 | 1 |
| 4i | 10.9908008588991 | 5.70902906222806 |
| 1/2+i | 0.384199141100893 | -0.351886205085204 |
| 5i | 17.1731262919990 | 9.24178645257085 |
| 1/2+5i/2 | 4.29312629199899 | 1.88178645257085 |

These are conditional-model coordinates, not measured A/E/C/W amplitudes.
Finite widths differ, and no continuum limit follows from exact feasibility.

## Three exact conditional-E4 identities

The displayed decimals obey the exact identities

`g(2i)=11/4`,
`g(4i)+g(1/2+i)=91/8`,
`g(5i)-g(1/2+5i/2)=322/25`.

For completeness, the prime coset identity is

`p^4 E4(p tau) + sum_(b=0)^(p-1) E4((tau+b)/p) = p(1+p^3) E4(tau)`.

It follows coefficient-by-coefficient from the q-series and
`sigma3(p n)+p^3 sigma3(n/p)=(1+p^3)sigma3(n)`, with the last term zero
when p does not divide n. The divisor identity follows immediately by splitting
off the power of p. Absolute convergence of the weight-four lattice sum also
gives the required S/T changes of basis.

At tau=i, p=2, the two inverse-related terms are 16 E4(2i) each,
while E4((1+i)/2)=-4 E4(i). Hence E4(2i)=11 E4(i)/16.
At tau=2i the same coset identity gives
`16 E4(4i)+E4(i)+E4(1/2+i)=18 E4(2i)`, yielding 91/8.
At tau=i, p=5, S/T reduction gives
`1250 E4(5i) - (625/2) E4(1/2+5i/2) -14 E4(i)=630 E4(i)`.
After area normalization this is `50(g(5i)-g(1/2+5i/2))=644`.
The exact statements do not depend on treating floating residuals as proofs.

## Minimal mechanism prediction coordinates

The two archived N50 four-vectors need not share a ray. Denote them x_i and
x_2i. A conditional **vector affine-E4** shape model has the frozen predictor

`x_j(N,tau) = c_j(N) [(1-t) x_j(50,i) + t x_j(50,2i)]`,

with separate vectors multiplying 1 and g. It neither assumes nor rehabilitates
a common scalar transport across A/E/C/W. The new-area factors c_j(N) are
additional assumptions or calibration parameters; they cannot be learned from
the N50 data alone. No exponent or value for them is selected here.

At N100, the ratio to a new 2i bridge cell cancels a coordinate's unknown factor:

`x_j(100,tau)/x_j(100,2i) = t+(1-t) x_j(50,i)/x_j(50,2i)`.

Source uncertainty and potentially unbounded Fieller intervals must be retained.
The N100 shear cell has a short extrapolation coefficient t=-0.351886, whereas
4i has t=5.709 and is more elongated. Thus **2i plus 1/2+i is the smallest
new-area bridge-and-discriminator pair**; all three N100 shapes add the exact
conditional sum constraint. This is a proposed prediction coordinate, not a
request to launch or a claim that its shape hypothesis survives production.
Replacing E4(tau) by E4(i Im(tau)) is a distinct height-only adversary, not a
legitimate modular-basis change for the shear cell.

## Scientific card and reproduction

- Changed space: N50 third-tau proposals are impossible; N100 is the minimal
  three-modulus fixed-map replacement, with a nondegenerate shear alternative.
- Not proved: E4 amplitudes, scale transport, field identity, or model survival.
- Observer/source/geometry: conditional P4(A,E,C,W), same rational map,
  square-site integer periods, fixed-area Smith pairs.
- Dependency: exact arithmetic; no new random block and no independent reuse
  of the N50 production estimates.
- Next useful output: evaluate surviving vector/clock hypotheses and their
  covariance-aware conditional N100 shear/2i predictions on the existing data.

```
python3 scripts/p267_third_geometry_certificate.py --json results/p267-third-geometry-feasibility/certificate.json
python3 -m unittest discover -s tests -p test_p267_third_geometry_certificate.py -v
```

Four focused tests passed. No Monte Carlo, remote connection, PR or comment.
