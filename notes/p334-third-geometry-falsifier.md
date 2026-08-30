# The paired balance times survive; Gaussian reinforcement does not

The third geometry was selected before any boundary flux was read. The rule
scanned HNF quotients in increasing order and required an honest face, no
physical quarter-turn symmetry, exactly two projective stabilizer orbits, and
distinct constant non-real `chi4` on both. The first candidate is

```text
P = [[7,2],[0,1]],  N=7,
ell_0=(0,1),        chi_0=(-7+24i)/25,
ell_1=(1,-3),       chi_1=( 7+24i)/25.
```

This lattice is not Gaussian-similar: any square-similarity lattice is
invariant under physical quarter turn, while the exact stabilizer test fails
here. The scan examined 29 geometry-only rows and never used birth, exit or
root information.

The full exact calculation has only 128 subsets and 448 directed additions.
The orbit-resolved continuity equation passes coefficientwise. The two net
balance times are

```text
p_1 = 4/7 = 0.571428571428571...,
p_0 = 0.592783237894885157837638571036...
```

Their separation `0.0213546665` is below the frozen N13/N17 envelope
`0.0409499231`: the close paired-zero timing prediction survives. At the
inherited `p_ref`, the upper orbit is only `3.72e-5` below its balance point;
its net is merely `1.09e-4` of source-plus-sink activity. The two root slopes
are `-6.67298` and `-1.67930`.

But the character geometry has changed decisively:

```text
Re(chi_0 conjugate(chi_1)) = 527/625 > 0.
```

The character contributions therefore reinforce when the scalar orbit nets
have the same sign and cancel when their signs differ. Exact phase topology:

```text
(0, 4/7)          reinforce
(4/7, p_0)        cancel
(p_0, 1)          reinforce.
```

This is the opposite of the Gaussian N13/N17 topology. Hence the frozen claim
“reinforcement only between the paired zeros” is falsified, even though the
paired timing zeros themselves survive.

The minimal replacement mechanism separates two layers:

1. source/sink balance dynamics generates a close pair of scalar net zeros;
2. the exact character Gram controls whether the between-zero interval
   reinforces or cancels.

This separation is an exact two-orbit theorem, not a fit. For real scalar
orbit-net currents `J_1,J_2` and fixed complex characters `chi_1,chi_2`,

```text
Re[(chi_1 J_1) conjugate(chi_2 J_2)]
  = Re[chi_1 conjugate(chi_2)] J_1 J_2
  = Gram(chi_1,chi_2) J_1 J_2.
```

The paired zeros control only the sign of `J_1 J_2`. Gaussian N13/N17 have
opposite-character `Gram<0`, whereas this HNF has `Gram=527/625>0`; the exact
Gram sign therefore flips the reinforcement/cancellation topology of the same
between-zero timing window.

This is also the clean bridge to #337: finite-field or ambient character
geometry supplies the Gram factor, while the #334 projective source/sink
process supplies `J_1 J_2`. They are separable factors, so a twist-character
change can flip alignment without changing the timing-zero mechanism.

The next falsifier is therefore parameter-free: for another two-orbit
quotient, positive character Gram must give between-zero cancellation and
negative Gram must give between-zero reinforcement. Failure rejects the
two-scalar-current reduction rather than merely changing a fitted share.
