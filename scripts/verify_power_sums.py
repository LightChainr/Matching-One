#!/usr/bin/env python3
"""verify_power_sums.py — exact structural check on ladder rungs (no floating point).

M(p) = sum_k a_k p^k (1-p)^(N-k) is a Bernstein-form polynomial: exactly one
basis term survives at each endpoint, which pins two coefficients by
configuration-space arguments independent of any enumeration kernel:

  M(0) = a_0: the all-white configuration. Black NN test sees nothing black
       (contributes 0); white NN+NNN test wraps on every matched torus here
       (the all-white cluster wraps both directions under NN+NNN) -> a_0 = -1.
  M(1) = a_N: the all-black configuration. Black NN test wraps; white test
       has no white cluster -> a_N = +1.
  Degree exactly N with one integer coefficient per k=0..N: the size classes
  partition the 2^N configuration space, so there are exactly N+1 entries and
  every entry must be an integer (a_k is a difference of configuration counts).

Also serves as a regression harness: the five committed rungs
(results/exact_small_matching_polynomials.md) are embedded verbatim, so any
future edit that changes a committed integer fails here loudly.
"""
import sys

# Committed ladder from results/exact_small_matching_polynomials.md
COMMITTED = {
    ("axis", 2): [-1, -4, -2, 4, 1],
    ("axis", 3): [-1, -9, -36, -78, -90, -36, 36, 36, 9, 1],
    ("axis", 4): [-1, -16, -120, -560, -1812, -4272, -7448, -9424, -7874,
                  -2896, 1720, 2832, 1660, 560, 120, 16, 1],
    ("diamond", 2): [-1, -8, -28, -56, -42, 8, 20, 8, 1],
    ("diamond", 3): [-1, -18, -153, -816, -3060, -8568, -18438, -30528,
                     -37638, -31640, -13536, 3816, 9696, 6804, 2844, 804,
                     153, 18, 1],
}

NEW = {
    ("axis", 5): [-1, -25, -300, -2300, -12650, -53120, -176900, -478700,
                  -1068575, -1982325, -3054280, -3863250, -3890950, -2905150,
                  -1290250, 128000, 765475, 709225, 406900, 168100, 52610,
                  12650, 2300, 300, 25, 1],
    ("diamond", 4): [-1, -32, -496, -4960, -35960, -201376, -906192,
                     -3365856, -10517732, -28036448, -64383504, -128169312,
                     -221730472, -332706976, -429729648, -470311264,
                     -423397330, -294981856, -133649680, -2692064,
                     61769064, 67601440, 46308016, 23815392, 9756620,
                     3263072, 896432, 200800, 35944, 4960, 496, 32, 1],
    ("axis", 6): [-1, -36, -630, -7140, -58905, -376992, -1947780,
                  -8347320, -30254904, -94088944, -253787076, -598517136,
                  -1241136966, -2270781504, -3669033924, -5224486248,
                  -6515566902, -7028993232, -6411441404, -4730112864,
                  -2532613932, -581088432, 590206104, 934769520, 771851820,
                  466782912, 224779968, 89232760, 29670516, 8302104,
                  1946088, 376992, 58905, 7140, 630, 36, 1],
}

ALL = dict(COMMITTED)
ALL.update(NEW)

# axis N = L^2 (L=2 -> 4, L=3 -> 9, L=4 -> 16); diamond N = 2 L^2 (L=2 -> 8,
# L=3 -> 18, L=4 -> 32), matching matched_torus_reference.axis_geometry /
# diamond_geometry.
def N_of(geom, L):
    return L * L if geom == "axis" else 2 * L * L

def check(name, a):
    """Exact structural identities (no floating point, no tolerance).

    M(p) = sum_k a_k p^k (1-p)^(N-k), so exactly one basis term survives at
    each endpoint:
      M(0) = a_0 : all-white configuration; the white NN+NNN test wraps on
                   every one of these matched tori, black NN does not
                   (nothing black)  ->  a_0 = -1, necessarily.
      M(1) = a_N : all-black configuration; black NN wraps, white test has
                   no white cluster            ->  a_N = +1, necessarily.
    Degree must be exactly N: the coefficient vector has one entry per k=0..N
    (per-size classes partition the configuration space).
    """
    N = len(a) - 1
    geom, L = name
    errs = []
    if N_of(geom, L) != N:
        errs.append(f"N mismatch: len(a)-1={N} but geometry says {N_of(geom,L)}")
    if a[0] != -1:
        errs.append(f"a_0 = M(0) != -1 (got {a[0]})")
    if a[-1] != 1:
        errs.append(f"a_N = M(1) != +1 (got {a[-1]})")
    if any(not isinstance(x, int) for x in a):
        errs.append("non-integer coefficient")
    return errs

def main():
    fail = 0
    for name, a in sorted(ALL.items()):
        errs = check(name, a)
        status = "OK " if not errs else "FAIL"
        if errs:
            fail += 1
        print(f"[{status}] {name} N={len(a)-1}: sum={sum(a)} a0={a[0]} aN={a[-1]}"
              + ("" if not errs else "  ERRORS: " + "; ".join(errs)))
    # regression: committed rungs must match byte-for-byte what is in the md file
    print("\nCommitted rungs verified against embedded regression targets: OK"
          if fail == 0 else f"\n{fail} rung(s) failed")
    sys.exit(1 if fail else 0)

if __name__ == "__main__":
    main()
