# Geometry-aware crosswalk of the N65/N85/N145 natural current

Status: no new simulation.  N65 and N85 select the model direction; both N145
components are held out and enter only the diagnostic score.

## Exact geometry dictionary

Every quotient is a Gaussian ideal square with exact modulus `tau=i`.
Consequently `E4(i)` is common to all six geometries and has no discriminating
rank.  For representation `a+ib`, the varying exact H4 covector is

```text
z_axis = ((a+ib)/sqrt(N))^4,
c4 = Re z_axis = cos(4 theta).
```

The charged A source also supplies a principled second descriptor.  On the
four lines of `P1(F3)`,

```text
q_A^2 = (1,1,0,0) = (u + H_F3)/2.
```

Thus the reflection-even charged response has the physical H4 covector and one
projective scalar `u/2=1/2`.  `Im z_axis` is retained in the certificate but is
not put into a reflection-even `K_A` model.

| N | representation | exact c4 | K_A |
|---:|---|---:|---:|
| 65 | 8+i | 3713/4225 | -0.04156368 |
| 65 | 7+4i | -2047/4225 | +0.02791698 |
| 85 | 9+2i | 4633/7225 | -0.01058760 |
| 85 | 7+6i | -6887/7225 | +0.01500502 |
| 145 | 12+i | 19873/21025 | -0.01072378 |
| 145 | 9+8i | -20447/21025 | +0.00958769 |

The machine result retains the full paired 2x2 covariance at each N and a
block-diagonal 6x6 covariance across the independent scale archives.

## Frozen model comparison

All models use the project area power `N^-13/8`; no exponent is fitted.

1. pure N law: one scalar radial column;
2. one H4 covector: `N^-13/8 c4`;
3. H4 plus charged/projective scalar:
   `N^-13/8 (beta_H c4 + beta_A/2)`.

| model | N65/N85 chi2 / df | N145 predictive chi2 / 2 | predicted N145 pair | pair residual |
|---|---:|---:|---:|---:|
| pure N | 23.134 / 3 | 42.593 | 0 | 6.476 SE |
| H4 geometry | 2.575 / 3 | 1.519 | 0.0152412 | 1.119 SE |
| H4 + A scalar | 2.454 / 2 | 1.029 | 0.0156005 | 1.014 SE |

The geometry-free radial direction already fails on N65/N85.  One exact H4
covector predicts the held-out N145 vector
`(-0.007512,+0.007729)` against observed
`(-0.010724,+0.009588)`.  Adding the A scalar improves source chi-square by
only `0.121`; it is not identified as a necessary second direction.

## Rotation versus scale curvature

Anchoring directly on the N85 pair gives:

```text
radial-only H4 target at N145      = 0.0107448,
geometry-aware H4 target at N145  = 0.0129231,
observed N145 pair                = 0.0203115.
```

The exact H4 pair covector grows by `1.20273` between the chosen N85 and N145
pairs.  It explains 22.8% of the central rebound above the radial-only target;
77.2% remains as a central scale-curvature term.  These percentages are an
accounting identity, not model probabilities.

With covariance, the geometry-aware N85-anchored remainder is only `1.453 SE`,
and the joint N65/N85 H4 direction scores `1.519/2` on N145.  The useful answer
is therefore two-level:

- geometry rotation is necessary and removes the apparent contradiction;
- central values still lean toward scale curvature, but an extra curvature
  parameter is not required at current precision.

This also explains why another arbitrary N would be a poor next experiment.

## Selected same-lineage geometry

The next geometry should be the N170 exact angle-flip child of N85:

```text
parent N85:  9+2i, 7+6i
multiply both by 1+i
child N170: 11+7i, 13+i
```

After canonical D4 representatives, the reflection-even H4 covector flips
sign exactly for both orientations, while the charged scalar remains `1/2`.
The area ratio is two and the fixed radial factor is `2^-13/8=0.3242099`.
N170 is also close enough to N145 to make a geometry comparison more direct.

No N170 archive currently stores the required projective sufficient
statistics.  The minimum missing row is the existing sparse schema:
`orientation,batch,samples,tau1,tau2,kind,ell_x,ell_y,count`.  This crosswalk
selects the geometry but does not authorize or run it.

## Reproduction

```bash
python3 scripts/crosswalk_natural_current_geometry.py \
  --json results/p337-natural-current-geometry-crosswalk/latest.json \
  --markdown results/p337-natural-current-geometry-crosswalk/latest.md

python3 -m unittest discover -s tests \
  -p 'test_crosswalk_natural_current_geometry.py'
```
