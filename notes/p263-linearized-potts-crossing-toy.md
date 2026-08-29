# Minimal two-channel tangent-crossing algebra

This is a bounded exact-first realization of Issue #263.  It does not attempt a Potts bootstrap.  It asks what crossing-vector shapes are forced by one simple-pole collision before any numerical conformal blocks are fitted.

## 1. Two opposite-residue channels

Write `epsilon=Q-1` and

```text
A_+ =  R/epsilon + a_+ + b_+ epsilon + ...,
A_- = -R/epsilon + a_- + b_- epsilon + ...,
Delta_+ = Delta_* + v_+ epsilon + w_+ epsilon^2/2 + ...,
Delta_- = Delta_* + v_- epsilon + w_- epsilon^2/2 + ... .
```

Expanding the two blocks gives the finite `Q=1` row

```text
G_0 = (a_+ + a_-) B
    + R(v_+-v_-) partial_Delta B,
```

and its `Q` tangent

```text
G_1 = (b_+ + b_-) B
    + [a_+v_+ + a_-v_- + R(w_+-w_-)/2] partial_Delta B
    + R(v_+^2-v_-^2)/2 partial_Delta^2 B.
```

This bookkeeping has one useful invariant: normalized OPE tangents and finite projector counterterms can move the `B` and `partial_Delta B` coefficients, but the displayed `partial_Delta^2 B` coefficient is fixed by the simple residue and the two dimension velocities.

There is also a hard regularity gate.  If the two dimensions do not coincide at `Q=1`, the term

```text
R [B(Delta_+(1))-B(Delta_-(1))]/epsilon
```

does not cancel.  Opposite Potts-projector residues do not turn arbitrary same-spin fields into a logarithmic pair.

## 2. Exact energy/two-cluster control

Use the #262 convention in which the singlet residue is `+R` and the `[2]` residue is `-R`, with `R=2J`.  The energy-family velocity follows from #261,

```text
v_energy = -9 sqrt(3)/(16 pi).
```

Together with the VJS positive-control gap

```text
d_Q(x_2cluster-x_energy)|_1 = sqrt(3)/pi,
```

this gives

```text
v_2cluster = 7 sqrt(3)/(16 pi),
R(v_energy-v_2cluster) = -R sqrt(3)/pi,
R(v_energy^2-v_2cluster^2)/2 = 3R/(16 pi^2).
```

The first number fixes the derivative-block coefficient already in the finite confluent correlator.  The second fixes the second-derivative-block coefficient in its `Q` tangent.  Neither depends on the normalized OPE tangent.

## 3. The spin-4 pair is an adversary, not this collision

The two current spin-4 candidates have

```text
Q4 epsilon:  x=21/4, x'=-9 sqrt(3)/(16 pi), singlet,
V_(2,2):     x=17/4, x'=-5 sqrt(3)/(16 pi), [2].
```

Their velocity gap is useful spectroscopy, but their dimensions differ by one at `Q=1`.  Therefore they are separate columns of the crossing matrix.  Applying the #262 opposite-residue formula directly to these two blocks would leave a `1/(Q-1)` divergence.  The `4:-6:3` Ward row from #252 remains the cheap primary-versus-descendant gate.

## 4. Where the other derivatives enter

For a four-point channel, the #250 normalized OPE tangent enters through

```text
kappa_p = d_Q log(C_12p C_34p).
```

Multiplication by `1+kappa_p epsilon+...` changes the finite Laurent coefficient `a_p`, hence the `B` and `partial_Delta B` rows above.  It does not manufacture the `partial_Delta^2 B` coefficient.

On the lattice, #258 supplies only

```text
M_i = Cov(G_i,T_measure).
```

The vector scored against crossing must be the reconstructed sum of measure score, finite confluent projector term, and explicit insertion derivative.  These pieces are recorded separately even though crossing constrains their sum.

## 5. A four-point rank score

At arbitrary sampled cross ratios, form the block-jet matrices

```text
M2 = [B, partial_Delta B],
M3 = [B, partial_Delta B, partial_Delta^2 B].
```

A rank-2 tangent obeys every three-point minor

```text
det[B, partial_Delta B, Y]=0.
```

A rank-3 second derivative can violate those minors but must obey the four-point relation

```text
det[B, partial_Delta B, partial_Delta^2 B, Y]=0.
```

The included stripped toy makes this explicit.  Put `t=(partial_Delta B)/B` and sample the symmetric coordinates

```text
t=(-2,-1,1,2),
z=logistic(t) approximately (0.119203,0.268941,0.731059,0.880797).
```

For `y=Y/B`, define

```text
r_left  = -2 y1 + 3 y2 - y3,
r_right =   - y2 + 3 y3 - 2 y4.
```

Then

```text
rank 2: r_left=r_right=0,
rank 3 with y=A+Bt+C t^2/2: r_left=r_right=-3C.
```

Thus the single frozen four-point null

```text
r_left-r_right=0
```

survives a genuine second-derivative block but rejects a generic fourth direction.  Two cross ratios cannot test rank two, three give its first null, and four distinguish a coherent rank-3 jet from unrestricted shape contamination.

The toy `t` coordinate is not asserted to be a Potts conformal block.  In the actual score, replace its columns by numerically evaluated `B`, `partial_Delta B`, and `partial_Delta^2 B` at the same declared cross ratios and propagate their full covariance.
