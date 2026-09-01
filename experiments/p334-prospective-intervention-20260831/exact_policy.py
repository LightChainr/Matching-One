"""Exact identities for the one declared affine intervention; no RNG or I/O."""
from fractions import Fraction


def policy(census, source, epsilon=Fraction(1, 8)):
    """census rows: (label, joint_safe_bool, e_first, e_second, L_first, L_second)."""
    if source not in (0, 1) or epsilon <= 0:
        raise ValueError("fixed physical source and positive epsilon required")
    d = len(census)
    if d == 0 or len({r[0] for r in census}) != d:
        raise ValueError("one complete vacant-label census required")
    classes = {}
    for row in census:
        if row[1]:
            classes.setdefault(tuple(row[2:4]), []).append(row)
    scores = {r[0]:Fraction(0) for r in census}
    for group in classes.values():
        n = len(group)
        marks = [r[4+source] for r in group]
        sm = sum(marks)
        for row in group:
            scores[row[0]] = Fraction(n*row[4+source]-sm, d)
    plus = {u:(1+epsilon*h)/d for u,h in scores.items()}
    minus = {u:(1-epsilon*h)/d for u,h in scores.items()}
    if min(*plus.values(), *minus.values()) <= 0:
        raise ValueError("positive-policy bound violated; never reduce epsilon after outcomes")
    assert sum(plus.values()) == sum(minus.values()) == 1
    for group in classes.values():
        assert sum(scores[r[0]] for r in group) == 0
        assert sum(plus[r[0]] for r in group) == sum(minus[r[0]] for r in group) == Fraction(len(group), d)
    weight = sum(max(h, 0) for h in scores.values())/d
    return scores, plus, minus, weight


def verify_finite_contrast(census, source, conditional_means, epsilon=Fraction(1, 8)):
    """Deterministic rational toy/old-census QA, not scientific sampling."""
    h, qp, qm, w = policy(census, source, epsilon)
    response = sum(h[u]*Fraction(conditional_means[u]) for u in h)/len(h)
    effect = sum((qp[u]-qm[u])*Fraction(conditional_means[u]) for u in h)
    assert effect == 2*epsilon*response
    if w:
        positive = sum(max(h[u], 0)*Fraction(conditional_means[u]) for u in h)/(len(h)*w)
        negative = sum(max(-h[u], 0)*Fraction(conditional_means[u]) for u in h)/(len(h)*w)
        assert w*(positive-negative) == response
    else:
        assert response == 0
    return response
