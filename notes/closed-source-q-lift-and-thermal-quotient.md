# The closed source selects a Q lift and an exact thermal-quotient tangent

**New consequence.** The already fixed source selects a particular positive
generic-Q family; its single-geometry matching root is exactly the
`Z_2-Q Z_0=0` critical-polynomial section of that family's unprojected
weights. It also fixes the Q tangent, including a rank counterterm which
cannot be erased by retuning density. A local occupied-edge perturbation
is the precise difference from the usual site random-cluster continuation.

This identifies the continuation rather than choosing a favorable Q lift
after observing a derivative. It does not prove Potts universality for
the resulting family. All statements through section 4 are finite-volume
identities; section 5 is an explicit conditional scaling prediction.

The earlier [Q-lift transport](https://github.com/LightChainr/Matching-One/blob/c34c838cd391fae890de72d732cc53279757263b/notes/generic-q-lift-semantics.md)
already distinguished the homology and critical-polynomial lifts. The new
step here is the fixed microscopic source selecting one complete measure,
path and sector normalization, not a new proof that generic lifts differ.

## 1. An exact three-parameter positive family

Use the same honest square torus, occupied site count K, occupied NN edge
count B, occupied component count C_B, and ambient image rank r. Write
`q=r-1`, `E=q^2`; neither observer is changed in the following family:

```text
P_(a,Q,eta)(A) = Z^-1 a^K exp(eta B) Q^(C_B-r/2),
a>0, Q>0, eta real.                                             (1)
```

Every weight is positive, including noninteger Q. The exact closed-source
law `Bernoulli_p exp(t S*)` is the one specified curve

```text
Q=exp(2t), eta=log Q, a=[p/(1-p)] Q^(-5/2),
S*=2C_B+2B-5K-r+2N+1.                                         (2)
```

The removed factor `(1-p)^N Q^(N+1/2)` is configuration-independent.
Equation (1) is an occupation-sum continuation, not a claim of a finite
number of colours at every real Q. The finite-colour construction at
integer `sqrt(Q)` in [the local gas note](closed-source-local-colour-gas.md)
provides a local realization of the same law with its topological factor.

The usual site random-cluster model has unprojected weights `a^K Q^C_B`.
In (1) it corresponds to eta=0 before applying the declared rank factor.
Our curve instead has `eta=log Q`: the extra `Q^B` is real, even though
both models reduce to the identical iid site law at Q=1, eta=0.

Wang et al. define the former model and report numerical agreement with
Potts exponents at Q=1.5,...,4, plus first-order signatures at larger Q.
Their result does not establish the phase diagram of the extra-edge
coupled curve (2). [Primary paper, section II.B](https://arxiv.org/abs/1411.4408v3).

## 2. The rank projection selects the critical-polynomial numerator

For each geometry separately, define the *unprojected*, unnormalized
sector sums of (1):

```text
A_j(a,Q,eta)=sum_(r=j) a^K exp(eta B) Q^C_B,       j=0,1,2.
```

Multiplying the projected normalizer by Q gives

```text
D=Q A_0+sqrt(Q) A_1+A_2,
<q>=(A_2-Q A_0)/D,
<E>=(A_2+Q A_0)/D.                                             (3)
```

Thus a single-geometry matching root is exactly `A_2=Q A_0`.
This is the weighted, rather than the unweighted `A_2=A_0`, homology
section. The selection follows from the pre-existing `Q^(-r/2)` factor;
no new observer coefficient has been optimized.

For the actual two-geometry observable the root remains

```text
0 = (1/2) sum_g (A_(2,g)-Q A_(0,g))/D_g.                         (4)
```

It is **not** a zero of the sum of unnormalized numerators. Each geometry
has its own partition function, and the pooling order in (4) is unchanged.

The source derivative of the numerator contains both the derivative of
the sector sums and the explicit term `-A_0 dQ`. In the projected measure
this same term is included automatically by `-r/2` in its log-Q score.
The two descriptions agree; dropping either counterterm changes the lift.

Deleting the rank factor instead gives `(A_2-A_0)/(A_0+A_1+A_2)` and the
single-geometry section `A_2=A_0`. This connects the two fixed source laws
directly to two different Q-lift sections. The overview team's
[opposite-sign tail result at fbbaa2aa](https://github.com/LightChainr/Matching-One/blob/fbbaa2aa/notes/topological-projection-reverses-global-u-tail.md)
already shows that this change has a real consequence for original U;
it is not repeated or counted as new numerical evidence here.

## 3. The root-adjusted response is an exact linear functional

At Q=1 use z=logit(p), denote the separately normalized pooled matching
mean by M(z), and keep Y(z), the original cos4 projection of E. Let
`M(z0)=0`, `M_z(z0)>0`, and

```text
U=A_N Y_z/M_z,       A_N=N^(13/8)/2.
```

For a fixed configuration statistic H, perturb by `exp(sH)` and define
`D_H U` as the total derivative along the *perturbed pooled root*.
All partial derivatives below retain the individual geometric normalizers:

```text
z_H = -M_H/M_z,
D_H U = A_N { (Y_zH+z_H Y_zz)/M_z
              -Y_z (M_zH+z_H M_zz)/M_z^2 }.                    (5)
```

Here e.g. `M_H` is the pooled covariance of q with H. Mixed derivatives
also include normalization, as derivatives of normalized expectations.
Equation (5) is linear in H. It annihilates constants and K exactly:

```text
D_1 U=0,       D_K U=0.                                        (6)
```

For K the perturbed law is simply z -> z+s; its root moves by -s and its
root-evaluated ratio is unchanged. This is a finite identity, not a claim
that every local source is equivalent to thermal retuning.

## 4. A non-arbitrary Q tangent and its explicit control direction

Combining (2), (5), and (6), at the common iid point gives

```text
(1/2) D_(S*) U = D_(C_B+B-r/2) U
               = D_(C_B-r/2) U + D_B U.                         (7)
```

The factor 1/2 is essential: `dQ/dt=2` at Q=1. The `-5K/2` part of the
fixed-p Q score cancels *only* after following the root as in (5).
The rank term does not cancel that way.

The first term on the last line is the root-adjusted Q tangent of the
rank-projected ordinary site random-cluster lift, at eta=0. The second
is the response to an independent local occupied-edge coupling at fixed
Q. Therefore the difference between two named Q continuations is exactly
one existing microscopic statistic, not an unspecified normalization
ambiguity or a new feature search.

At the iid root, B itself has an exact density-free representative. With
`T_2=sum_(ij NN) (n_i-p0)(n_j-p0)` and the same fixed p0 in all terms,
four-regularity gives `B=T_2+4p0 K-2N p0^2`. Hence `D_B U=D_(T_2)U`
by (6). No assumption about CFT scaling is needed for this local control.

For reference, the fixed-density rank score is already determined by
the existing q/E state, geometry by geometry:

```text
partial_s <q>_(exp(sr)) |0 = <E>-<q>^2,
partial_s <E>_(exp(sr)) |0 = <q>(1-<E>).                         (8)
```

Thermal differentiation of (8), followed by (5), gives the full rank
counterterm. The pooled root does not set each geometric `<q>` to zero.
There is no need to introduce a new topological observable to retain it.

## 5. A two-direction prediction that can fail

Suppose a smooth critical family through (1) lies on a single Potts branch,
and the same nonzero field dominates original U for a fixed shape pair:

```text
U_N(Q,eta)=A(Q,eta) N^[13/8-(x(Q)-2)/2] (1+small corrections).     (9)
```

The fixed-Q local coordinate eta can alter the amplitude and irrelevant
corrections, but in this *specified hypothesis* does not alter x. The
root and normalizations are those of (3)-(5). For equal area multiplier
c at the same shapes define

```text
L_S(N,c)=[(D_(S*)U/U)_(cN)-(D_(S*)U/U)_N]/log(c),
L_B(N,c)=[(D_B U/U)_(cN)-(D_B U/U)_N]/log(c).
```

If the corrections and their first parameter derivatives tend to zero,
equations (7)-(9) require the joint limits

```text
L_S -> -x'_Q(1),             L_B -> 0.                           (10)
```

This gives a necessary control for interpreting the source response as
field spectroscopy. A persistent nonzero L_B rules out the stated
single-field, fixed-Q-universal-exponent description, even if L_S happens
to resemble a favored Q velocity. It does not by itself identify a new
field. Mixtures, Jordan logarithms, vanishing leading amplitudes, and a
different critical family are specified alternatives to (9), not reasons
to add a descriptor until (10) holds.

The exact source identity (7) holds regardless of (9). No finite pair of
sizes is automatically an asymptotic exclusion; a production comparison
must state its correction allowance before reveal. This note does not
launch a block, prescribe adaptive sample increases, or reopen stopped
P154/P334 experiments.

The [companion weak-source calculation](closed-source-weak-colour-spectroscopy.md)
retains the later ordinary four-leg selection zero. A nonzero derivative
of a vanishing four-leg endpoint overlap, if the declared Q continuation
permits it, instead gives `R_N=alpha+gamma_T log N+beta sqrt(N)+...`, with
`gamma_T=9 sqrt(3)/(16 pi)`. After removing that known thermal logarithm,
resolved successive increments at equal area ratio c have ratio sqrt(c).
This is not the old constant four-leg velocity and does not resurrect a
regular overlap already excluded by the selection theorem. A source-
activated Jordan logarithm can also shift a constant slope; slope constancy
alone is therefore not a field identification.

## Scientific card

- **Mechanism changed:** the closed source now has a fully declared Q
  continuation and weighted homology root; its difference from the usual
  site-RC tangent is exactly the local B response.
- **Not established:** Potts universality, a critical line, single-field
  dominance, or a measured cross-size limit.
- **Observer/sector/source/geometry:** original pooled q/E and U; ordinary
  spin4 shape projection; fixed S*, with B as the specified Q-path control;
  honest equal-area square-site shape pairs.
- **Dependency:** algebraic consequence of the same action and rank
  projection, not an additional random-data evidence group.
- **Next discriminator:** a joint, predeclared cross-size comparison of
  L_S and L_B after the finite-size correction model has been specified.
  Exact failure of (10)'s second limit would stop the single-field
  interpretation, not trigger a search for a more favorable Q lift.

Literature was retrieved using the arXiv skill. Only the primary
[site-RC paper](https://arxiv.org/abs/1411.4408v3) is used here; the identities
and conditional control above are derived from the repository's fixed law.
