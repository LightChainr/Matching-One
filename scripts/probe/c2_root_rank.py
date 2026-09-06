#!/usr/bin/env python3
"""C2: root-count algebra of the #549 fork language.

Two parts.

Part 1 re-verifies the #549 closed form from its published protocol constants
(BASE_SAFE counts and the EXIT_A/EXIT_B successor tables of #435's two N16
gadgets).  The fork protocol has a single "root" (one shared update vertex)
and two conditionally independent clones; the success probability is affine in
the class coordinate a (the number of A-type gadgets in the bank of k).

Part 2 tests the mechanism claim with an abstract exact gadget model:
    - a bank of k gadgets, class a = number of A-type gadgets;
    - a test with d independent roots (d sequential shared updates, each on a
      fresh gadget) and one clone per root;
    - success = every root/clone event succeeds, per-gadget success
      probability h_tau depending only on the gadget type tau in {A,B}.

Claim: the success probability is then a polynomial of degree d in a (d = root
count = language depth), so d = 1 gives rank <= 2 (the #549 case), and d = k
gives a Vandermonde-style rank k+1 response matrix.  Rank growth needs *depth*
(more independent roots), not clone width.

All arithmetic exact (Fraction).  No sampling.
"""

from fractions import Fraction
from math import comb

# ---------------------------------------------------------------------------
# Part 1: #549 closed form and its successor-sum derivation
# ---------------------------------------------------------------------------

BASE_SAFE = (1, 7, 18, 20, 8, 0, 0, 0, 0)   # safe m-subset counts of the #435 N16 gadget
EXIT_A = (1, 1, 1, 2, 2, 3, 3)              # safe-successor exit counts, gadget A
EXIT_B = (1, 2, 2, 2, 2, 2, 2)              # safe-successor exit counts, gadget B


def fork_probability_closed(k: int, a: int) -> Fraction:
    return Fraction(343 * k**3 - 182 * k**2 + 25 * k + 4 * a, 8 * k * (8 * k - 1) ** 2)


def fork_probability_by_successor_sum(k: int, a: int) -> Fraction:
    """Direct protocol count: one root vertex, two clones, safe-second counts."""
    exits = [EXIT_A] * a + [EXIT_B] * (k - a)
    total = 0
    for local in exits:
        for x in local:
            s = 7 * k - x              # safe second vertices after a root of "exit x"
            total += s * s
    return Fraction(total, 8 * k * (8 * k - 1) ** 2)


def verify_part1(k_max: int = 12) -> None:
    for k in range(1, k_max + 1):
        for a in range(k + 1):
            assert fork_probability_closed(k, a) == fork_probability_by_successor_sum(k, a)
        for a in range(k):
            gap = fork_probability_closed(k, a + 1) - fork_probability_closed(k, a)
            assert gap == Fraction(1, 2 * k * (8 * k - 1) ** 2)
    print(f"part 1: #549 closed form == successor-sum protocol count for k=1..{k_max}  [PASS]")


# ---------------------------------------------------------------------------
# Part 2: root-count algebra (abstract exact gadget model)
# ---------------------------------------------------------------------------

def poly_degree_in_a(test_success: "callable", k: int) -> int:
    """Finite-difference degree of a symmetric class-probability in a.

    Returns the degree if the values on a = 0..k come from a polynomial of
    that degree (checked by (deg+1)-th finite differences vanishing).
    """
    vals = [test_success(k, a) for a in range(k + 1)]
    diffs = list(vals)
    degree = 0
    for order in range(1, k + 1):
        diffs = [diffs[i + 1] - diffs[i] for i in range(len(diffs) - 1)]
        if all(d == 0 for d in diffs):
            return order - 1
    return k


def make_droot_test(h_A: Fraction, h_B: Fraction, d: int):
    """Success probability of a d-root test in the abstract model.

    Protocol: choose d distinct gadgets uniformly (ordered, without
    replacement); on each, its clone must land in a good point, whose fraction
    is h_tau for a gadget of type tau.  The class is a = #A-type gadgets.
    """
    def prob(k: int, a: int) -> Fraction:
        if d > k:
            return Fraction(0)
        total = Fraction(0)
        # iterate over ordered distinct choices: equivalent to summing over
        # subsets weighted by per-gadget h, times ordering factor.
        # P = sum over ordered (g1..gd) distinct of prod h(gj) / (k)_d
        # For type counts: choose t A-positions in the ordered tuple and the
        # rest B.  Number of ordered tuples with exactly t A's =
        #   C(d,t) * (a)_t * (k-a)_{d-t}   (falling factorials)
        from math import factorial
        num = Fraction(0)
        for t in range(d + 1):
            if t <= a and d - t <= k - a:
                falling_a = factorial(a) // factorial(a - t)
                falling_b = factorial(k - a) // factorial(k - a - (d - t))
                num += Fraction(comb(d, t) * falling_a * falling_b, 1) \
                    * (h_A ** t) * (h_B ** (d - t))
        denom = factorial(k) // factorial(k - d)   # (k)_d ordered distinct
        return num / Fraction(denom, 1)
    return prob


def response_matrix(k: int, tests, a_range):
    rows = []
    for a in a_range:
        rows.append([tests[d](k, a) for d in range(len(tests))])
    return rows


def rank_fraction(rows) -> int:
    """Exact rank by fraction Gaussian elimination."""
    rows = [[Fraction(x) for x in r] for r in rows]
    m = len(rows)
    if m == 0:
        return 0
    ncols = len(rows[0])
    rank = 0
    for col in range(ncols):
        pivot = None
        for r in range(rank, m):
            if rows[r][col] != 0:
                pivot = r
                break
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        pv = rows[rank][col]
        rows[rank] = [x / pv for x in rows[rank]]
        for r in range(m):
            if r != rank and rows[r][col] != 0:
                f = rows[r][col]
                rows[r] = [rows[r][c] - f * rows[rank][c] for c in range(ncols)]
        rank += 1
        if rank == m:
            break
    return rank


def main() -> None:
    verify_part1()

    h_A, h_B = Fraction(2, 3), Fraction(1, 3)   # distinguishable types
    k = 5
    print(f"\npart 2: abstract d-root tests, k={k}, h_A={h_A}, h_B={h_B}")
    for d in range(1, k + 1):
        test = make_droot_test(h_A, h_B, d)
        deg = poly_degree_in_a(test, k)
        print(f"  d={d} roots: degree of P(a) = {deg}   (claim: degree = d)")
    print("\n  response matrix for tests d = 0..k (rows a = 0..k; d=0 = trivial test):")
    tests = [make_droot_test(h_A, h_B, d) for d in range(0, k + 1)]
    M = response_matrix(k, tests, a_range=range(k + 1))
    print(f"  exact rank(M) = {rank_fraction(M)}  (claim: rank = k+1 = {k + 1})")

    # two-type degeneracy control: h_A == h_B must collapse rank to 1
    tests0 = [make_droot_test(Fraction(1, 2), Fraction(1, 2), d) for d in range(1, k + 1)]
    M0 = response_matrix(k, tests0, a_range=range(k + 1))
    print(f"  degeneracy control (h_A=h_B): exact rank = {rank_fraction(M0)}  (claim: 1)")

    # #549 as the d=1 instance with a *single* root but two clones: still affine
    vals = [fork_probability_closed(k, a) for a in range(k + 1)]
    diffs = [vals[i + 1] - vals[i] for i in range(k)]
    print(f"\n  #549 single-root/two-clone fork at k={k}: 1st differences "
          f"{'all equal' if len(set(diffs)) == 1 else 'NOT affine'} -> affine in a "
          f"(degree 1, rank <= 2)")


if __name__ == "__main__":
    main()
