# Q→1 matching tangent as an exact intertwiner obstruction

There is an exact operator hiding between the two overly broad descriptions “S/D are only empirical covectors” and “S/D are physical parity fields.”

Take the doubled finite configuration space `H_G + H_Ghat` and define the complement/species exchange

```text
J |G,A>    = |Ghat,A^c>,
J |Ghat,B> = |G,B^c>.
```

Then `J^2=I` exactly. Let

```text
K = diag(k_G(A), k_Ghat(B))
```

be the generic-`Q` cluster tangent of the weight `exp((log Q)K)`. Its exact exchange components are

```text
S = (K+JKJ)/2,
D = (K-JKJ)/2,
J S J = S,
J D J = -D.
```

Thus S/D really are parity eigenoperators—but of the **doubled exchange superoperator** `Ad_J`, not automatically fields of either single lattice theory. More sharply,

```text
[K,J] = 2 D J,
exp(tK)J - J exp(tK) = 2 exp(tS) sinh(tD) J,  t=log Q.
```

So the matching-odd tangent `D` is exactly half the first-order failure of the bare complement interface to pull through the generic-`Q` cluster weight. At `Q=1` all weights equal one and `J` intertwines trivially; its first derivative is the obstruction. This realizes the “derivative defect” literally without pretending the defect remains topological away from `Q=1`.

## Positive control and no-go

For edge FK on `C3` and its planar dual, the raw tangent difference is the local Euler term `2-|A|`. Subtracting it makes `D=0` configuration by configuration. This is the exact positive control: the repaired interface genuinely intertwines at first order.

For the square-face site pair `C4 -> K4`, no counterterm depending only on occupancy `|A|` can do this. At `|A|=2`:

```text
A={0,1} adjacent:  k_C4(A)-k_K4(A^c)=0,
A={0,2} opposite:  k_C4(A)-k_K4(A^c)=1.
```

This is a finite exact no-go for the entire scalar local-density counterterm class. It does not rule out a richer connectivity-state interface, but it proves that the missing state must resolve geometry/topology inside a fixed occupancy sector.

The same witness separates this object from the #258 measure derivative. Both configurations have identical Bernoulli score—and score zero at `p=1/2`—while their matching `Q` tangent differs. Therefore the matching tangent cannot be absorbed into a critical-manifold/occupancy score; it belongs in the projector/defect or explicit-insertion terms of the three-part derivative ledger.

Finally, `Ad_J` parity is independent of the `S_Q` colour representation in #257. The matching `D` above is colour singlet. The exact `[2]` selection zero for `V_(2,+/-2)` therefore does not remove this channel, while singlet thermal `Q4` remains eligible. “Matching odd” and “Potts charged” are distinct quantum numbers.

## Boundary

This proves a finite doubled-space intertwiner identity and an occupancy-counterterm obstruction. It does not yet construct a bounded-local transfer seam for the infinite square-site model, identify empirical `P4[S']/P4[D']` with these operators, or imply a continuum Jordan block. The next local construction must reproduce `D` while carrying enough connectivity information to distinguish the adjacent/opposite witness.
