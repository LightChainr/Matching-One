#!/usr/bin/env python3
"""Exact thermal density-warp obstructions and one common-channel area.

Coefficients are rational, low-degree first. This computes algebraic
functionals, not statistical significance of empirical polynomial means.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import json
from math import comb
from pathlib import Path


def trim(p):
    p = list(map(F, p)) or [F(0)]
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def add(p, q):
    return trim([(p[i] if i < len(p) else 0)+(q[i] if i < len(q) else 0)
                 for i in range(max(len(p), len(q)))])


def scale(p, c):
    return trim([F(c)*x for x in p])


def sub(p, q):
    return add(p, scale(q, -1))


def mul(p, q):
    out = [F(0)]*(len(p)+len(q)-1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i+j] += a*b
    return trim(out)


def deriv(p):
    return trim([i*p[i] for i in range(1, len(p))])


def primitive(p):
    return [F(0)]+[F(a, i+1) for i, a in enumerate(p)]


def value(p, x):
    out = F(0)
    for a in reversed(p):
        out = out*F(x)+a
    return out


def integral(p, a=0, b=1):
    h = primitive(p)
    return value(h, b)-value(h, a)


def compose(p, phi):
    out = [F(0)]
    for a in reversed(p):
        out = add(mul(out, phi), [a])
    return out


def divmod_poly(p, q):
    r, q = trim(p), trim(q)
    if q == [0]:
        raise ZeroDivisionError
    out = [F(0)]*max(1, len(r)-len(q)+1)
    while r != [0] and len(r) >= len(q):
        degree, coefficient = len(r)-len(q), r[-1]/q[-1]
        out[degree] = coefficient
        r = sub(r, [F(0)]*degree+scale(q, coefficient))
    return trim(out), r


def gcd_poly(p, q):
    p, q = trim(p), trim(q)
    while q != [0]:
        _, r = divmod_poly(p, q)
        p, q = q, r
    return scale(p, 1/p[-1])


def roots_open_unit(p):
    """Distinct roots in (0,1); endpoints must be nonzeros."""
    p = trim(p)
    if value(p, 0) == 0 or value(p, 1) == 0:
        raise ValueError("remove endpoint factors first")
    seq = [p, deriv(p)]
    if seq[-1] == [0]:
        return 0
    while True:
        _, r = divmod_poly(seq[-2], seq[-1])
        if r == [0]:
            break
        seq.append(scale(r, -1))
    def changes(x):
        signs = [1 if value(q, x) > 0 else -1 for q in seq if value(q, x) != 0]
        return sum(a != b for a, b in zip(signs, signs[1:]))
    return changes(0)-changes(1)


def tangent_gate(d, r, alpha=1):
    """Endpoint-fixed regular R=alpha(vD)' iff the rational v has no poles
    on [0,1] and vanishes at both endpoints. A finite warp is a different gate.
    """
    if trim(d) == [0] or F(alpha) == 0:
        raise ValueError("The source polynomial and amplitude must be nonzero")
    h = primitive(r)
    denominator = scale(d, alpha)
    common = gcd_poly(h, denominator)
    numerator, _ = divmod_poly(h, common)
    denominator, _ = divmod_poly(denominator, common)
    endpoint_poles = [p for p in (0, 1) if value(denominator, p) == 0]
    interior_denominator = denominator
    for p in (0, 1):
        while value(interior_denominator, p) == 0:
            interior_denominator, rem = divmod_poly(interior_denominator, [-F(p), F(1)])
            assert rem == [0]
    poles = roots_open_unit(interior_denominator)
    endpoints = ([value(numerator, p)/value(denominator, p) for p in (0, 1)]
                 if not endpoint_poles else None)
    return {"mass_residual": integral(r), "velocity_numerator": numerator,
            "velocity_denominator": denominator, "interior_pole_count": poles,
            "endpoint_poles": endpoint_poles, "velocity_endpoints": endpoints,
            "regular_endpoint_fixed_tangent": poles == 0 and endpoints == [0, 0],
            "boundary": "Only a linear-generator statement; cannot reject a finite density warp from this alone."}


def density_pullback(d, phi):
    return mul(deriv(phi), compose(d, phi))


def omega(a, e):
    """Oriented area of the cumulative A/E path; reparametrization invariant."""
    return integral(sub(mul(primitive(a), e), mul(primitive(e), a)))


def clock_moments(a, e, order=6):
    """J_m=integral (F_A/M_A)^m E dp, m=0,...,order.

    Algebraically defined when M_A != 0. Interpretation as moments in a
    monotone cumulative clock additionally requires A to have one sign.
    """
    mass = integral(a)
    if mass == 0:
        raise ValueError("A has zero mass; no normalized cumulative clock")
    clock, power = scale(primitive(a), 1/mass), [F(1)]
    out = []
    for _ in range(order+1):
        out.append(integral(mul(power, e)))
        power = mul(power, clock)
    return out


def bernstein_to_power(b):
    n = len(b)-1
    return trim([comb(n, k)*sum((-1)**(k-i)*comb(k, i)*F(b[i])
                                 for i in range(k+1)) for k in range(n+1)])


def encode(x):
    if isinstance(x, F):
        return str(x)
    if isinstance(x, dict):
        return {k: encode(v) for k, v in x.items()}
    if isinstance(x, (tuple, list)):
        return [encode(v) for v in x]
    return x


def supplied_area_gate(data):
    read = bernstein_to_power if data.get("basis", "power") == "bernstein" else trim
    da, de, ua, ue = [read(data[k]) for k in ("D_A", "D_E", "U_A", "U_E")]
    ra = F(data["r_A"] if "r_A" in data else data["alpha"])
    re = F(data["r_E"] if "r_E" in data else data["alpha"])
    jd, ju = clock_moments(da, de), clock_moments(ua, ue)
    return encode({"source_omega": omega(da, de), "target_omega": omega(ua, ue),
                   "common_density_warp_area_null": omega(ua, ue)-ra*re*omega(da, de),
                   "mass_null_A": integral(ua)-ra*integral(da),
                   "mass_null_E": integral(ue)-re*integral(de),
                   "clock_moment_nulls_m0_to_m6": [u-re*d for u, d in zip(ju, jd)],
                   "boundary": "Exact functionals of supplied polynomials, not uncertainty-calibrated rejections. A must be single-sign for a monotone clock interpretation. Full same-batch covariance, including fitted amplitudes, must be propagated externally."})


def certificate():
    p, one_minus_p = [F(0), F(1)], [F(1), F(-1)]
    envelope = mul(p, one_minus_p)
    da = mul(envelope, [F(-1, 3), F(1)])
    de = mul(envelope, [F(1), F(1)])
    v = scale(envelope, F(1, 5))
    good_r = deriv(mul(v, da))
    bad_r = mul(envelope, [F(-1, 2), F(1)])
    assert integral(bad_r) == 0
    bad_flux = value(primitive(bad_r), F(1, 3))
    assert bad_flux == F(-1, 81)
    good = tangent_gate(da, good_r)
    bad = tangent_gate(da, bad_r)
    assert good["regular_endpoint_fixed_tangent"] and not bad["regular_endpoint_fixed_tangent"]
    assert bad["interior_pole_count"] == 1
    nozero = tangent_gate(envelope, bad_r)
    assert nozero["regular_endpoint_fixed_tangent"]

    # A multiple zero requires vanishing order, not just H(z)=0.
    double = mul(envelope, mul([F(-1, 2), F(1)], [F(-1, 2), F(1)]))
    h = mul(mul(envelope, envelope), [F(-1, 2), F(1)])
    multiple = tangent_gate(double, deriv(h))
    assert value(h, F(1, 2)) == 0 and multiple["interior_pole_count"] == 1

    phi = add(p, scale(envelope, F(1, 4)))
    phi_e = sub(p, scale(envelope, F(1, 4)))
    assert value(phi, 0) == 0 and value(phi, 1) == 1
    assert roots_open_unit(deriv(phi)) == 0 and value(deriv(phi), F(1, 2)) > 0
    ua, ue = density_pullback(da, phi), density_pullback(de, phi)
    area_equal = omega(ua, ue)-omega(da, de)
    assert area_equal == 0
    finite_h = primitive(sub(ua, da))
    finite_flux = value(finite_h, F(1, 3))
    finite_tangent = tangent_gate(da, sub(ua, da))
    assert finite_flux != 0 and not finite_tangent["regular_endpoint_fixed_tangent"]
    separate_e = density_pullback(de, phi_e)
    separate_gap = omega(ua, separate_e)-omega(da, de)
    assert integral(ua) == integral(da) and integral(separate_e) == integral(de)
    assert separate_gap != 0

    # A nonvanishing clock and independent, including negative, amplitudes.
    clock_a, clock_e = [F(1), F(1)], [F(-1, 3), F(1)]
    ra, re = F(-2, 7), F(3, 5)
    clock_ua = scale(density_pullback(clock_a, phi), ra)
    clock_ue = scale(density_pullback(clock_e, phi), re)
    jd, ju = clock_moments(clock_a, clock_e), clock_moments(clock_ua, clock_ue)
    moment_gaps = [u-re*d for u, d in zip(ju, jd)]
    assert moment_gaps == [0]*7
    independent_amplitude_gap = omega(clock_ua, clock_ue)-ra*re*omega(clock_a, clock_e)
    assert independent_amplitude_gap == 0

    simple_d, simple_g = [F(-1, 3), F(1)], [F(-5, 6), F(2)]
    source_lobes = [integral(simple_d, 0, F(1, 3)), integral(simple_d, F(1, 3), 1)]
    target_lobes = [integral(simple_g, 0, F(5, 12)), integral(simple_g, F(5, 12), 1)]
    tv_gap = sum(map(abs, target_lobes))-sum(map(abs, source_lobes))
    assert integral(simple_d) == integral(simple_g) and tv_gap == F(17, 72)
    return encode({
        "schema": "matching-one/p267-thermal-warp-invariants/v1",
        "status": "exact theoretical certificates; no production fit or simulation",
        "scalar_vs_density": {"scalar": "G=D composed with phi; R=v D'",
                              "density": "G=phi' (D composed with phi); R=(vD)'"},
        "global_tangent_criterion": "H=integral_0^p R; v=H/(alpha D). The reduced denominator must have no roots on [0,1], and v(0)=v(1)=0.",
        "regular_example": good,
        "zero_area_but_singular_example": {"D": da, "R": bad_r, "zero": F(1, 3),
                                           "flux_at_zero": bad_flux, "simple_pole_residue": F(-1, 18), "gate": bad},
        "same_residual_zero_free_source": nozero,
        "multiple_zero_counterexample": {"D": double, "H": h, "gate": multiple},
        "finite_warp_not_linear_generator": {"phi": phi, "D": da, "G": ua,
                                              "old_zero_flux": finite_flux, "tangent_gate": finite_tangent,
                                              "finite_density_warp_verified": True},
        "finite_TV_obstruction": {"D": simple_d, "G": simple_g,
                                  "source_signed_lobe_areas": source_lobes, "target_signed_lobe_areas": target_lobes,
                                  "same_total_mass": integral(simple_d), "TV_gap": tv_gap},
        "common_channel_area": {"definition": "Omega(D)=integral(F_A D_E-F_E D_A)dp",
                                "finite_null": "Omega(U)-r_A*r_E Omega(D)=0; equal amplitudes give alpha^2",
                                "same_phi_gap": area_equal, "different_individually_valid_phi_gap": separate_gap,
                                "source_omega": omega(da, de)},
        "single_sign_clock": {"D_A": clock_a, "D_E": clock_e, "phi": phi,
                              "r_A": ra, "r_E": re,
                              "independent_amplitude_area_null": independent_amplitude_gap,
                              "clock_moment_nulls_m0_to_m6": moment_gaps,
                              "relation": "When mass nulls vanish, J1(U)-r_E*J1(D)=(Omega(U)-r_A*r_E*Omega(D))/(2*r_A*M_A)",
                              "completeness": "All moments determine the signed E measure in q=F_A/M_A; finitely many moment nulls are necessary, not sufficient."},
        "first_variation_area": "delta Omega=2 integral(H_A D_E-H_E D_A)+M_A H_E(1)-M_E H_A(1)",
        "statistics_boundary": "Exact nonzero polynomials from empirical means are not statistical no-go certificates; use the same-stream covariance or batch jackknife for TV, lobes and Omega.",
    })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, help="D_A,D_E,U_A,U_E, r_A/r_E (or common alpha); basis power or bernstein")
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    result = supplied_area_gate(json.loads(args.input.read_text())) if args.input else certificate()
    text = json.dumps(result, indent=2)+"\n"
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(text)
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
