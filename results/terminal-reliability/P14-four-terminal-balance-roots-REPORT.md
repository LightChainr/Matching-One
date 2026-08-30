# P14 four-terminal necessary balance-root screen

## Result

All 27 frozen connected one-internal-vertex candidates have exactly one root
in `(0,1)` of

```text
B_G(p) = P_G(0123;p) - P_G(0|1|2|3;p).
```

The roots were obtained from exact integer power-basis coefficients and
isolated with an exact-rational Sturm sequence to intervals narrower than
`2^-96`.

The root nearest the descriptive square-site value is

```text
graph 4:5:0001001101
edges (0,4), (1,4), (2,3), (3,4)
B(p) = -1 + p + 3p^2 - 5p^3 + 3p^4
p_balance = 0.6062392719967533...
distance = 0.0134932212071533...
```

It is a tree with four bridges, two articulation vertices and nontransitive
terminals. The fully terminal-transitive four-spoke star is the next closest
at `0.6245191881991872...`, but its hub is an articulation and all four edges
are bridges. A numerically attractive root in this family is therefore a
reducibility effect, not a credible four-terminal replacement cell.

Exactly one candidate passes all preregistered structural gates:

```text
graph 4:5:0111111011
edges (0,2), (0,3), (0,4), (1,2), (1,3), (1,4), (2,4), (3,4)
structure: four-terminal wheel W5 (cyclic order 0,2,1,3 plus hub 4)
terminal automorphism order: 8; terminal-transitive
bridges/articulations: 0/0
zero crossing pair for the compatible cyclic order: 01|23
p_balance = 0.2979305190327643...
distance = 0.2948155317568357...
```

Its exact balance polynomial is

```text
-1 + 4p - 24p^3 + 96p^4 - 176p^5 + 164p^6 - 76p^7 + 14p^8.
```

Thus the bounded homogeneous family has a sharp split: roots near the square-
site number occur only in structurally reducible/asymmetric graphs, while the
unique structurally credible disk candidate is the dense wheel and lives in a
very different parameter range.

## Mechanism decision

Do not enlarge the same homogeneous one-internal census merely to chase a
closer scalar root. The useful continuation is the wheel's two natural edge
orbits: rim and spoke probabilities. A planar dual exchanges those orbits, so
the next exact object should be the full two-parameter partition vector on the
duality line, not another one-parameter `P_all=P_none` ranking. Only after a
periodic tiling/comparison map is specified can that object enter Issue #14's
rigorous-bound route.

## Boundary

This is a necessary scalar balance screen. It supplies no embedding
certificate, full partition-duality equality, periodic construction,
stochastic domination, comparison baseline, critical point or rigorous bound.

## Reproduction

```bash
python3 scripts/p14_four_terminal_balance_roots.py \
  --output results/terminal-reliability/p14-four-terminal-balance-roots.json
uv run --with pytest python -m pytest -q \
  tests/test_p14_four_terminal_balance_roots.py
```

Focused plus inherited exact-corpus result: `27 passed`.
