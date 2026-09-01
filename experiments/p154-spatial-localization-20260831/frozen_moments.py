"""Unchanged binomial_moments extracted from Matching-One 764595ea5c838c110e416382a3a90e2ecf7297bb."""
import math
import numpy as np
from scipy.special import gammaln

def binomial_moments(sums, samples, p, n):
    """Every permutation supplies all K; derivatives integrate Binomial weights."""
    k = np.arange(n + 1, dtype=float)
    log_b = gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)
    log_b += k * math.log(p) + (n - k) * math.log1p(-p)
    weights = np.exp(log_b)
    score = k / p - (n - k) / (1 - p)
    wp = weights * score
    wpp = weights * (score * score - k / p**2 - (n - k) / (1 - p)**2)
    z, zp, zpp = weights.sum(), wp.sum(), wpp.sum()
    mean = weights @ sums / (samples * z)
    first = (wp @ sums / samples - mean * zp) / z
    second = (wpp @ sums / samples - mean * zpp - 2 * first * zp) / z
    return mean, first, second, {'binomial_mass': float(z), 'permutations': samples}
