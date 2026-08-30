# Adjacent annihilator formal asymptotics

This bounded Issue #47 oracle isolates the exponent mapping used by the
Mertens–Ziff adjacent-size criterion.  With

```text
M_L(pc) = a L^(-13/4) [1+c L^(-q)+...],
M'_L(pc) = b L^(3/4) [1+...],
```

linearizing `L^(13/4)M_L=(L-1)^(13/4)M_(L-1)` gives

```text
delta p = (a c/b) L^(-(4+q)) [q/4+O(1/L)].
```

The script obtains the full checked prefix of the bracket from the exact formal
quotient

```text
- [1-(1-x)^(-q)] / [1-(1-x)^4],  x=1/L,
```

after stripping the common leading `x`.  Fraction generalized-binomial
coefficients cover both integer and half-integer candidates.  The frozen map is
`q={3/2,2,3,4,6}` to `w={11/2,6,7,8,10}` with leading coefficient `q/4`.

Run:

```text
python3 scripts/adjacent_annihilator_asymptotics.py
python3 -m unittest tests/test_adjacent_annihilator_asymptotics.py -v
```

## Boundary

The certificate assumes the displayed single-correction expansion and a
nonzero leading slope amplitude.  It reads no historical or production roots,
does not select an exponent, estimate `p_c`, or validate the asymptotic regime.
The held-out numerical challenge in Issue #47 remains open.
