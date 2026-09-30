#!/usr/bin/env python3
"""One small-k exact control of parity encoding, tensor rank and readout TV.

Standard library only. k=2,3,4; no physical graph enumeration, old-census
rerun, Monte Carlo, cloud access, or full future-query rank claim.
Print JSON to stdout; the geometric and all-k proofs are in the note.
"""

from collections import Counter
from fractions import Fraction
from itertools import product
import json
from math import comb


LOGICAL = ((0, 0), (1, 0), (0, 1), (1, 1))
EPSILONS = (Fraction(1, 10), Fraction(1, 4))


def qstr(value):
    value = Fraction(value)
    return f"{value.numerator}/{value.denominator}"


def exact_rank(rows):
    basis = []
    for source in rows:
        row = [Fraction(x) for x in source]
        for pivot, vector in basis:
            factor = row[pivot]
            if factor:
                row = [x - factor*y for x, y in zip(row, vector)]
        pivot = next((i for i, value in enumerate(row) if value), None)
        if pivot is not None:
            leading = row[pivot]
            basis.append((pivot, [x/leading for x in row]))
    return len(basis)


def monotone_truth_tables(k):
    # A monotone function on k variables consists uniquely of monotone
    # lower/upper sections f0 <= f1 on k-1 variables. Avoid enumerating
    # all 2^(2^k) Boolean tables or all encoders.
    tables = [0, 1]
    for dimension in range(1, k+1):
        shift = 1 << (dimension-1)
        tables = [low | (high << shift)
                  for low in tables for high in tables if low & ~high == 0]
    return tables


def parity_difference(table, k):
    signed = sum((-1)**x.bit_count() for x in range(1 << k)
                 if table >> x & 1)
    return Fraction(signed, 1 << (k-1))


def analyze(k):
    full = (1 << k)-1
    classes = [[x for x in range(1 << k) if x.bit_count() % 2 == bit]
               for bit in (0, 1)]
    code_supports = [set(product(classes[u], classes[v])) for u, v in LOGICAL]
    assert all(len(support) == 1 << (2*k-2) for support in code_supports)
    assert sum(map(len, code_supports)) == len(set.union(*code_supports)) == 4**k

    marginal_checks = 0
    for support in code_supports:
        for subset in range(full):  # every proper subset of physical blocks
            counts = Counter((x & subset, y & subset) for x, y in support)
            assert len(counts) == 4**subset.bit_count()
            assert set(counts.values()) == {4**(k-subset.bit_count()-1)}
            marginal_checks += 1

    queries = list(product(range(3), repeat=k))
    rows = []
    moment_checks = 0
    for (u, v), support in zip(LOGICAL, code_supports):
        row = []
        for query in queries:
            a = sum(1 << i for i, q in enumerate(query) if q == 1)
            b = sum(1 << i for i, q in enumerate(query) if q == 2)
            actual = Fraction(sum(x & a == a and y & b == b for x, y in support),
                              len(support))
            predicted = Fraction(1, 1 << (a.bit_count()+b.bit_count()))
            if a == full:
                predicted *= 1 + (-1)**(k+u)
            if b == full:
                predicted *= 1 + (-1)**(k+v)
            assert actual == predicted
            row.append(actual)
            moment_checks += 1
        rows.append(row)
    rank = exact_rank(rows)
    assert rank == 3
    assert all(rows[0][j]-rows[1][j]-rows[2][j]+rows[3][j] == 0
               for j in range(len(queries)))
    selected = [queries.index((q,)*k) for q in range(3)]
    product_tv = max(abs(rows[0][j]-rows[1][j]) for j in range(len(queries)))
    assert product_tv == Fraction(1, 1 << (k-1))

    tables = monotone_truth_tables(k)
    differences = [parity_difference(table, k) for table in tables]
    bound = Fraction(comb(k-1, (k-1)//2), 1 << (k-1))
    assert max(map(abs, differences)) == bound
    threshold_results = []
    for threshold in range(1, k+1):
        table = sum(1 << x for x in range(1 << k)
                    if x.bit_count() >= threshold)
        delta = parity_difference(table, k)
        assert delta == Fraction((-1)**threshold * comb(k-1, threshold-1),
                                 1 << (k-1))
        if abs(delta) == bound:
            threshold_results.append(threshold)

    noise_results = []
    for epsilon in EPSILONS:
        attenuation = (1-2*epsilon)**k
        # Direct finite BSC summation, not the Fourier formula under test.
        noisy = []
        for parity in (0, 1):
            distribution = []
            for y in range(1 << k):
                mass = sum(epsilon**((x ^ y).bit_count()) *
                           (1-epsilon)**(k-(x ^ y).bit_count())
                           for x in classes[parity]) / len(classes[parity])
                formula = Fraction(1, 1 << k) * (
                    1 + (-1)**(parity+y.bit_count())*attenuation)
                assert mass == formula
                distribution.append(mass)
            assert sum(distribution) == 1
            noisy.append(distribution)
        full_tv = sum(abs(a-b) for a, b in zip(*noisy))/2
        assert full_tv == attenuation
        parity_error = sum(noisy[0][y] for y in range(1 << k)
                           if y.bit_count() % 2)
        assert parity_error == (1-attenuation)/2
        noisy_tvs = []
        for table, delta in zip(tables, differences):
            noisy_delta = sum(noisy[0][y]-noisy[1][y]
                              for y in range(1 << k) if table >> y & 1)
            assert noisy_delta == attenuation*delta
            noisy_tvs.append(abs(noisy_delta))
        assert max(noisy_tvs) == attenuation*bound
        noise_results.append({
            "epsilon": qstr(epsilon),
            "all_noisy_shares_TV": qstr(full_tv),
            "parity_decoder_error": qstr(parity_error),
            "max_single_monotone_TV": qstr(max(noisy_tvs)),
        })

    return {
        "k": k,
        "parity_class_size": len(classes[0]),
        "joint_support_size_per_logical_state": len(code_supports[0]),
        "four_supports_pairwise_disjoint": True,
        "proper_block_marginal_checks": marginal_checks,
        "tensor_columns": len(queries),
        "exact_moment_checks": moment_checks,
        "exact_tensor_response_rank": rank,
        "alternating_row_relation": [1, -1, -1, 1],
        "constant_allU_allV_columns": [[qstr(row[j]) for j in selected] for row in rows],
        "max_single_product_query_TV_for_u": qstr(product_tv),
        "monotone_boolean_functions_checked": len(tables),
        "max_single_monotone_TV": qstr(bound),
        "best_single_monotone_equal_prior_error": qstr((1-bound)/2),
        "maximizing_monotone_functions": sum(abs(delta) == bound for delta in differences),
        "maximizing_weight_thresholds": threshold_results,
        "hypothetical_independent_share_BSC": noise_results,
    }


def main():
    cases = [analyze(k) for k in (2, 3, 4)]
    result = {
        "schema": "matching-one.stochastic-parity-block-encoding.v1",
        "contract": {
            "logical_state_order": [list(state) for state in LOGICAL],
            "encoding": "independent uniform U,V subject to xor(U)=u and xor(V)=v",
            "encoder_randomness": "private; no side channel revealing the sampled shares",
            "response_space": "full tensor of local (1,U_i,V_i)",
            "general_physical_future_query_rank_upper_bound_proved_here": False,
            "monotone_TV_scope": "one Boolean output; not a transcript or adaptive controller",
            "BSC_scope": "independent errors on observed shares, not derived physical instruction errors",
        },
        "arithmetic": "exact integers and fractions.Fraction",
        "cases": cases,
        "total_tensor_moments_checked": sum(case["exact_moment_checks"] for case in cases),
        "total_monotone_functions_checked": sum(case["monotone_boolean_functions_checked"] for case in cases),
        "failures": 0,
        "not_run": ["old deterministic census", "physical graph enumeration", "Monte Carlo",
                    "cloud jobs", "whole-machine audit"],
    }
    print(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
