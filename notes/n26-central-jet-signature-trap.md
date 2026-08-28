# N=26 central-jet trap: a globally false Beta law can mimic kappa3/kappa5

**Status:** post-target exact analysis of the open PR #152 N=26 polynomial. This note does not alter any frozen target.

The exact N=26 self-matching result supplies a useful warning for the universal-scaling-function program.

PR #152 gives the exact power-basis matching polynomial

```text
M(p) = -1 + 156 p^5 - 338 p^6 + 260 p^7 - 260 p^8 - 338 p^9
       + 1144 p^10 + 3536 p^11 - 13702 p^12 + 15628 p^13
       - 3016 p^14 - 10088 p^15 + 11492 p^16 - 5798 p^17
       + 1482 p^18 - 156 p^19.
```

At `p=1/2`, exact differentiation gives

```text
M'     = 342927 / 65536
M'''   = -1629693 / 8192
M^(5)  = 2441985 / 128
```

and therefore

```text
kappa3(N=26)
 = M'''/(M')^3
 = -159532265242624 / 114894152000433
 = -1.3885151025095062...

kappa5(N=26)
 = M^(5)/(M')^5
 = 591381954649765670111150080
   /121602861284634003091801113
 = 4.863224009717393...
```

Now compare the completely explicit `Beta(7,7)` centered law. The reliability-signature note gives

```text
kappa3[Beta(7,7)]
 = -4194304/3006003
 = -1.3953093193852435...

kappa5[Beta(7,7)]
 = 43980465111040/9036054036009
 = 4.867220241908278...
```

The relative differences are only

```text
kappa3: about 0.487%
kappa5: about 0.0821%.
```

Yet the same Beta(7,7) law is **exactly impossible globally** on N=26: it predicts zero successful configurations at occupation `k=5`, while the exact geometry already has 78 successful configurations in the canonical `F=(1+M)/2` direction.

So N=26 provides an exact example of the following phenomenon:

> A wrong global threshold-distribution family can reproduce low-order metric-free central derivatives extremely well.

This matters directly for the continuum program around `kappa3` and `kappa5`.

## Consequences

1. A near-simple value of `kappa3` is not sufficient to identify a scaling-function class.
2. Even adding `kappa5` can leave a strongly wrong global family almost indistinguishable at accessible precision.
3. Quantiles, tails, multi-u values, or the full activation signature provide genuinely independent identification power.
4. The Gaussian/majority limit `(-pi/2, 3pi^2/4)` should be treated as a shape null, but it should be challenged jointly with global profile observables rather than only derivative invariants.
5. Issue #122 (full standardized threshold distribution) and the multi-u program #119 become more valuable in light of this exact finite counterexample.

The key point is methodological, not that Beta(7,7) is a serious N=26 model after its frozen failure. The exact example demonstrates how central-jet agreement can coexist with a categorical tail/support mismatch.
