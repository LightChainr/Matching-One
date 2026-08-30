# Exact two-orbit HNF atlas through index 12

The atlas applies one flux-blind gate to every HNF quotient of index 4 through
12: honest face, connected primal graph, no physical quarter-turn symmetry,
exactly two projective line orbits under the lattice D4 stabilizer, and
distinct nonzero effective orbit characters. Six of 119 HNFs pass:

```text
N7: [[7,2],[0,1]], [[7,3],[0,1]], [[7,4],[0,1]], [[7,5],[0,1]]
N8: [[8,3],[0,1]], [[8,5],[0,1]]
```

The complete reveal is only 1,024 subset states and 3,840 directed boundary
edges. All orbit continuity identities pass coefficientwise. The N8 reflection
orbits contain two lines; their line-resolved birth/exit coefficient tables
are exactly equal, so replacing the two characters by their orbit average is
an exact compression rather than a modelling assumption.

## Timing answer

Every included HNF has exactly one simple net-flux root in each orbit. There is
no root-count counterexample through N12. But closeness under the frozen
N13/N17 envelope is not universal:

```text
N7 separation = 0.0213546664663   pass
N8 separation = 0.0506069555622   fail (> 0.0409499231130)
```

Thus the minimal close-pair counterexample is `[[8,3],[0,1]]`. The robust
bounded statement is one balance time per orbit, not a universal narrow
separation.

## Character-Gram answer

The six HNFs reduce to two exact mechanism signatures after symmetry copies:

```text
N7: Gram =  527/625, reinforce / cancel / reinforce
N8: Gram =     -7/25, cancel / reinforce / cancel
```

For every coefficient of every full curve,

```text
Re[(chi1 J1) conjugate(chi2 J2)]
  = Gram(chi1,chi2) J1 J2.
```

All six polynomial residuals are identically zero. Gram sign therefore gives
a complete stratification of cooperation topology; this part is exact and
does not depend on root closeness.

The #337/#334 separation is now atlas-level rather than a single example:
character geometry supplies the Gram sign, while projective source/sink timing
supplies `J1 J2`. Beyond N12, root multiplicity and spacing remain dynamical
questions, but any quotient passing the same two-orbit compression gate must
obey the Gram-sign topology on every interval.
