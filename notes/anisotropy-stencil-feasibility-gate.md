# An exact stencil gate for the improved-action search

Status: exact Phase 0 geometry result for Issue 106. The parent control program remains open.

## Question

The proposed improved-action program seeks an exactly-critical percolation family whose leading
spin-4 amplitude can be tuned through zero. The cheapest candidate is the inhomogeneous square-bond
family with horizontal and vertical probabilities

```text
p_h=t,
p_v=1-t.
```

The equality `p_h+p_v=1` is its exact critical surface. Before simulating that family, this note asks
a narrower question: can its edge stencil make the microscopic fourth angular moment proposed in the
issue vanish?

For an edge vector `v=(dx,dy)`, define

```text
A4(v) = |v|^4 exp(4 i theta) = (dx+i dy)^4.
```

For nonnegative local weights `w_e`, the stencil proxy is

```text
A4 = sum_e w_e (dx_e+i dy_e)^4.
```

Integer vectors make this an exact integer/rational calculation; no trigonometric or floating-point
rounding is needed.

## Proposition 1: the axis-only family cannot cross zero

The four unit axial edges have directions `+/-x,+/-y`. Their fourth moments are all `+1`:

```text
(+/-1)^4 = 1,
(+/-i)^4 = 1.
```

If the horizontal pair has per-edge weight `w_h>=0` and the vertical pair has per-edge weight
`w_v>=0`, then

```text
A4_axis = 2 w_h + 2 w_v.
```

This is real and strictly positive unless all edge weights vanish. Therefore no nontrivial
horizontal/vertical-only stencil with nonnegative weights can tune this proxy through zero.

Two common proxy choices make the obstruction explicit on the exact critical surface:

```text
probability weights:       A4(t)=2t+2(1-t)=2,
Bernoulli-variance weights: A4(t)=4t(1-t).
```

The second expression vanishes only at `t=0,1`, where the random two-dimensional model degenerates.
The conclusion is not tied to either proxy: the general nonnegative formula already proves the gate.

There is a second experimental problem with `t!=1/2`: unequal horizontal and vertical probabilities
break C4 to C2 and introduce lower-symmetry response. Thus this family is neither a zero-stencil nor a
clean C4-preserving interpolation.

## Proposition 2: cancelling the positive axis shell needs a negative spin-4 phase

Assume reflection-complete edge orbits, so imaginary fourth moments cancel. If the positive axial
shell is retained and every added orbit has `cos(4 theta)>=0`, a nonnegative weighted sum cannot
vanish. Hence cancelling that axial contribution requires an orbit in a sector with

```text
cos(4 theta) < 0.
```

The integer diagonal is the minimal square-lattice example:

```text
(1+i)^4 = -4,
(1-i)^4 = -4.
```

For C4-complete unit-axis and integer-diagonal shells with per-edge weights `w_axis,w_diagonal`,

```text
A4_axis     =  +4 w_axis,
A4_diagonal = -16 w_diagonal,
A4_total    = 4 w_axis - 16 w_diagonal.
```

The exact geometric cancellation is therefore

```text
w_axis / w_diagonal = 4.
```

The factor four depends on the declared `|v|^4` and integer-step normalization: a diagonal has squared
length two and hence fourth radial weight four. Normalizing all directions to unit length would change
the numerical ratio. The manifest freezes the integer-step convention so this cannot drift silently.

## Decision

Do not spend simulation time scanning the exactly-critical axes-only square-bond family for a zero of
this microscopic proxy. Retain the larger Issue 106 program, but require its next candidate to have:

1. an exact criticality argument (isoradial or star--triangle transport);
2. at least one negative-phase oblique edge orbit;
3. a derived, not guessed, mapping from edge probabilities/geometry to the physical spin-4 coupling;
4. only then, a frozen H4 amplitude test.

The axis/diagonal `4:1` result is a design certificate, not a claim that a corresponding independent
bond model is exactly critical.

## Evidence boundary

What is exact here:

- the inhomogeneous square-bond critical surface used as input;
- every integer fourth moment;
- the axes-only nonnegative-weight no-go;
- the necessary negative-phase condition for cancelling a retained positive axis shell under
  reflection symmetry;
- the axis/diagonal cancellation ratio for the declared normalization.

What is not established:

- that probability, variance, or any other convenient proxy weight equals the renormalized spin-4
  coupling;
- that an exactly-critical mixed axis/diagonal family exists at the geometric ratio;
- that the measured torus H4 amplitude vanishes when the microscopic proxy vanishes.

The exact critical surface and the broader isoradial/star--triangle route are anchored in
Grimmett--Manolescu, `arXiv:1105.5535` and `arXiv:1204.0505`.
