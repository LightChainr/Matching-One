#!/usr/bin/env python3
"""Exact Bell4 equality-pattern algebra; no colour or lattice enumeration."""
import json
import sympy as s

Q = s.Symbol("Q")


def partitions(prefix=(0,)):
    if len(prefix) == 4:
        yield prefix
    else:
        for value in range(max(prefix) + 2):
            yield from partitions(prefix + (value,))


patterns = tuple(partitions())


def kernels(p):
    a, b, c, d = p
    gate = int(a != b and c != d)
    p0 = gate / (Q * (Q - 1))
    p2 = s.Rational(1, 2) * gate * (
        int(a == c and b == d) + int(a == d and b == c)
        - sum((a == c, a == d, b == c, b == d)) / (Q - 2)
        + 2 / ((Q - 1) * (Q - 2))
    )
    return p0, p2


raw = {p: kernels(p) for p in patterns}
averaged = {}
for p in patterns:
    x = kernels(p)
    y = kernels(p[1:] + p[:1])
    averaged[p] = tuple(s.cancel((a + b) / 2) for a, b in zip(x, y))


def pairing(left, right):
    return s.factor(sum(s.ff(Q, max(p) + 1) * left(p) * right(p)
                        for p in patterns))


def compatible(external, actual):
    return all(external[i] != external[j] or actual[i] == actual[j]
               for i in range(4) for j in range(i))


def label(p):
    names = "NESW"
    return "|".join("".join(names[i] for i, x in enumerate(p) if x == b)
                    for b in range(max(p) + 1))


inner = {
    "K2bar_squared": pairing(lambda p: averaged[p][1], lambda p: averaged[p][1]),
    "K0bar_squared": pairing(lambda p: averaged[p][0], lambda p: averaged[p][0]),
    "K2bar_K0bar": pairing(lambda p: averaged[p][1], lambda p: averaged[p][0]),
    "Kreg_squared": pairing(lambda p: sum(averaged[p]), lambda p: sum(averaged[p])),
    "K2_quarter_turn_overlap": pairing(lambda p: raw[p][1],
                                        lambda p: kernels(p[1:] + p[:1])[1]),
}
out = {"exact_pairings": {k: str(v) for k, v in inner.items()},
       "K2bar_Q4": str(inner["K2bar_squared"].subs(Q, 4)),
       "K2bar_Q1_residue": str(s.limit((Q - 1) * inner["K2bar_squared"], Q, 1)),
       "Kreg_Q1_value": str(inner["Kreg_squared"].subs(Q, 1)),
       "Kreg_Q1_derivative": str(s.diff(inner["Kreg_squared"], Q).subs(Q, 1)),
       "closures": []}
for external in patterns:
    closure = sum(s.ff(Q, max(p) + 1) * sum(averaged[p])
                  for p in patterns if compatible(external, p))
    beta = s.factor(closure / Q ** (max(external) + 1))
    out["closures"].append({"partition": label(external), "beta": str(beta),
                            "Q1_value": str(beta.subs(Q, 1)),
                            "Q1_derivative": str(s.diff(beta, Q).subs(Q, 1))})
print(json.dumps(out, indent=2))
