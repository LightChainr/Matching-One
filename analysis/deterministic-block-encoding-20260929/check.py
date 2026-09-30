#!/usr/bin/env python3
"""One exact census of four labelled states encoded into two corner blocks.

Print the deterministic result JSON. This checks only the stated full tensor
query algebra, not reachability of those columns by a physical rank experiment.
No dependencies, simulation, floating-point rank, or whole-machine enumeration.
"""

from collections import Counter
from itertools import product
import json
from math import comb, factorial, gcd


CORNERS = ((0, 0), (1, 0), (0, 1), (1, 1))


def integer_rank(rows):
    """Gaussian elimination over Q using nonzero integer row operations.

    Eliminate with b[pivot]*row - row[pivot]*b. Primitive normalization
    controls coefficient growth; no reduction modulo a prime is used.
    """
    basis = []
    for source in rows:
        row = list(source)
        for pivot, b in basis:
            if row[pivot]:
                factor, leading = row[pivot], b[pivot]
                row = [leading * x - factor * y for x, y in zip(row, b)]
                divisor = 0
                for x in row:
                    divisor = gcd(divisor, x)
                if divisor:
                    row = [x // divisor for x in row]
        pivot = next((i for i, x in enumerate(row) if x), None)
        if pivot is not None:
            if row[pivot] < 0:
                row = [-x for x in row]
            basis.append((pivot, row))
    return len(basis)


def tensor_row(word):
    first, second = CORNERS[word % 4], CORNERS[word // 4]
    return tuple(x * y for x in (1, *first) for y in (1, *second))


def main():
    rows = tuple(tensor_row(word) for word in range(16))
    rank_counts = Counter()
    image_rank_counts = Counter()
    injective_cardinality_counts = Counter()
    deficient_informative_block_counts = Counter()
    examples = {
        (0, 1, 2, 3): "one_informative_block",
        (0, 5, 10, 15): "direct_copy",
        (0, 1, 4, 5): "split_bits",
        (0, 5, 6, 3): "parity_second_block",
    }
    example_results = {}
    checked = 0
    for code in product(range(16), repeat=4):
        actual = integer_rank(rows[word] for word in code)
        image_size = len(set(code))
        sizes = (len({word % 4 for word in code}),
                 len({word // 4 for word in code}))
        exception = image_size == 4 and sizes in ((4, 1), (1, 4))
        predicted = 3 if exception else image_size
        if actual != predicted:
            raise AssertionError((code, sizes, actual, predicted))
        rank_counts[actual] += 1
        image_rank_counts[image_size, actual] += 1
        if image_size == 4:
            injective_cardinality_counts[sizes, actual] += 1
            if exception:
                deficient_informative_block_counts[sizes.index(4)] += 1
        if code in examples:
            example_results[examples[code]] = {
                "corner_word_ids": list(code), "exact_rank": actual,
                "local_image_cardinalities": list(sizes),
            }
        checked += 1

    # Independent combinatorial count: choose the distinct words, then a
    # surjection from the four labelled states; deficient injections have
    # a unique informative block, a constant other corner, and 4! labels.
    surjections = {1: 1, 2: 14, 3: 36, 4: 24}
    distinct_word_counts = {m: comb(16, m) * n for m, n in surjections.items()}
    deficient_injections = 2 * 4 * factorial(4)
    expected_ranks = {
        1: distinct_word_counts[1],
        2: distinct_word_counts[2],
        3: distinct_word_counts[3] + deficient_injections,
        4: distinct_word_counts[4] - deficient_injections,
    }
    assert checked == 16 ** 4
    assert dict(rank_counts) == expected_ranks
    assert sum(deficient_informative_block_counts.values()) == deficient_injections

    result = {
        "schema": "matching-one.deterministic-block-encoding.v1",
        "contract": {
            "logical_states": 4,
            "physical_blocks": 2,
            "corner_order": [list(corner) for corner in CORNERS],
            "word_id": "first_corner + 4*second_corner",
            "response": "(1,u0,v0) tensor (1,u1,v1)",
            "columns": ["1", "u1", "v1", "u0", "u0*u1", "u0*v1",
                        "v0", "v0*u1", "v0*v1"],
            "query_language": "full algebraic tensor span",
            "physical_product_query_reachability_verified_here": False,
        },
        "arithmetic": "exact integer Gaussian elimination over Q",
        "assignments_checked": checked,
        "rank_counts": {str(k): v for k, v in sorted(rank_counts.items())},
        "by_distinct_codewords_and_rank": [
            {"distinct_codewords": m, "rank": r, "assignments": n}
            for (m, r), n in sorted(image_rank_counts.items())
        ],
        "injective_assignments": distinct_word_counts[4],
        "injective_rank3": deficient_injections,
        "injective_rank4": expected_ranks[4],
        "injective_by_local_cardinality": [
            {"local_image_cardinalities": list(sizes), "rank": r,
             "assignments": n}
            for (sizes, r), n in sorted(injective_cardinality_counts.items())
        ],
        "rank3_informative_block_counts": {
            str(k): v for k, v in sorted(deficient_informative_block_counts.items())
        },
        "unlabelled_injective_codebooks": {
            "total": comb(16, 4), "rank3": deficient_injections // factorial(4),
            "rank4": expected_ranks[4] // factorial(4),
        },
        "classification_mismatches": 0,
        "combinatorial_count_crosscheck": "passed",
        "examples": example_results,
        "general_k_rank3_injective_count": "24*k*4^(k-1), k>=1",
        "scope": "finite algebraic census supporting the all-k proof in the note",
        "not_claimed": [
            "all tensor columns are physically reachable in every source contract",
            "a stochastic-encoder obstruction",
            "a full-machine positive versus linear dimension theorem",
            "optimal noisy error probability",
            "novelty certification",
        ],
    }
    print(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
