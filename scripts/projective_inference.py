#!/usr/bin/env python3
"""Covariance-weighted tests of fixed signed-real model subspaces.

Full-rank inference retains the original GLS/Fieller calculation. Rank deficiency
is not automatically missing information or deterministic information: the
caller must say which. Default mode rejects singular/truncated covariance.
``covariance_mode='exact_support'`` imposes its nullspace as a known affine
support, profiles the remaining amplitudes, and computes the constrained df.
Empirical low rank from too few batches does NOT license exact_support.

These are Gaussian-reference/asymptotic tests, not exact finite-sample coverage
claims for an estimated covariance. A signed line is not a nonnegative ray.
"""
from __future__ import annotations
from typing import Any, Sequence
from mpmath import mp

DEFAULT_DPS = 50
DEFAULT_RELATIVE_TOLERANCE = '1e-12'


def _as_matrix(rows):
    n = len(rows)
    if not n or any(len(row) != n for row in rows):
        raise ValueError('covariance must be nonempty and square')
    result = mp.matrix([[mp.mpf(str(x)) for x in row] for row in rows])
    if any(not mp.isfinite(x) for x in result):
        raise ValueError('covariance entries must be finite')
    return result


def _as_vector(values):
    result = mp.matrix([mp.mpf(str(x)) for x in values])
    if any(not mp.isfinite(x) for x in result):
        raise ValueError('vector entries must be finite')
    return result


def _spectral(covariance, relative_tolerance):
    tol = mp.mpf(relative_tolerance)
    if not 0 < tol < 1:
        raise ValueError('relative_tolerance must lie in (0,1)')
    matrix = _as_matrix(covariance)
    n = matrix.rows
    for i in range(n):
        for j in range(i+1,n):
            if abs(matrix[i,j]-matrix[j,i]) > mp.mpf('1e-30')*(
                    abs(matrix[i,j])+abs(matrix[j,i])+1):
                raise ValueError('covariance is not symmetric')
            matrix[j,i] = matrix[i,j]
    values, vectors = mp.eigsy(matrix)
    eigenvalues = [values[i] for i in range(n)]
    scale = max(abs(v) for v in eigenvalues)
    # PSD validation is NOT the user-selected statistical rank cutoff.
    roundoff = 100*mp.eps*scale
    if min(eigenvalues) < -roundoff:
        raise ValueError('covariance is not positive semidefinite')
    largest = max(eigenvalues)
    kept = [i for i,v in enumerate(eigenvalues) if v > tol*largest]
    pinv = mp.zeros(n,n)
    for k in kept:
        inverse = 1/eigenvalues[k]
        for i in range(n):
            for j in range(n):
                pinv[i,j] += inverse*vectors[i,k]*vectors[j,k]
    condition = largest/min(eigenvalues[k] for k in kept) if kept else None
    return pinv, kept, condition, eigenvalues, vectors, roundoff


def spectral_pseudo_inverse(covariance, relative_tolerance=DEFAULT_RELATIVE_TOLERANCE):
    """Algebraic pseudoinverse and diagnostics; not by itself a test."""
    pinv,kept,condition,values,_,_ = _spectral(covariance,relative_tolerance)
    return pinv,len(kept),condition,values


def _least_squares_and_kernel(matrix, target, tol):
    """Minimum-norm fit, null basis and identifiable rank via SVD."""
    width = matrix.cols
    if matrix.rows == 0:
        return mp.zeros(width,1), mp.eye(width), 0
    if width == 0:
        return mp.zeros(0,1), mp.zeros(0,0), 0
    U,s,Vt = mp.svd(matrix,full_matrices=True)
    largest = max(s) if len(s) else mp.mpf(0)
    retained = [i for i in range(len(s)) if s[i] > tol*largest]
    solution = mp.zeros(width,1)
    for i in retained:
        coefficient = mp.fsum(U[j,i]*target[j] for j in range(matrix.rows))/s[i]
        for j in range(width):
            solution[j] += Vt[i,j]*coefficient
    free = [i for i in range(width) if i not in retained]
    null = mp.matrix(width,len(free))
    for col,i in enumerate(free):
        for j in range(width):
            null[j,col] = Vt[i,j]
    return solution,null,len(retained)


def _exact_support_fit(y,directions,kept,values,vectors,tol,roundoff):
    n,width = y.rows,len(directions)
    dropped = [i for i in range(n) if i not in kept]
    if any(abs(values[i]) > roundoff for i in dropped):
        raise ValueError('a numerical cutoff is not an exact support: lower relative_tolerance')
    V = mp.matrix([[directions[j][i] for j in range(width)] for i in range(n)])
    U0 = mp.matrix([[vectors[i,j] for i in range(n)] for j in dropped])
    constraints, rhs = U0*V,U0*y
    a0,Z,constraint_rank = _least_squares_and_kernel(constraints,rhs,tol)
    support_tol = mp.power(10,-mp.dps//2)
    if mp.norm(constraints*a0-rhs) > support_tol*max(1,mp.norm(y),mp.norm(V*a0)):
        raise ValueError('model is incompatible with exact covariance support')
    W = mp.matrix(len(kept),n)
    for row,k in enumerate(kept):
        for i in range(n):
            W[row,i] = vectors[i,k]/mp.sqrt(values[k])
    if Z.cols and kept:
        fit,_,fitted_rank = _least_squares_and_kernel(W*V*Z,W*(y-V*a0),tol)
        amplitudes = a0+Z*fit
    else:
        fitted_rank, amplitudes = 0,a0
    residual = W*(y-V*amplitudes)
    statistic = mp.fsum(x*x for x in residual)
    degrees = len(kept)-fitted_rank
    if degrees == 0:
        # Saturation is an algebraic zero, not a chi-square fluctuation.
        scale = max(1, mp.norm(W*y), mp.norm(W*V*amplitudes))
        if mp.norm(residual) > support_tol*scale:
            raise ArithmeticError('saturated support solve lost numerical accuracy')
        statistic = mp.mpf(0)
    return statistic,degrees,amplitudes,constraint_rank,fitted_rank,support_tol


def subspace_residual(observed: Sequence[Any], covariance: Sequence[Sequence[Any]],
                      basis: Sequence[Sequence[Any]],
                      relative_tolerance=DEFAULT_RELATIVE_TOLERANCE,
                      dps=DEFAULT_DPS, *, covariance_mode='strict') -> dict[str,Any]:
    if covariance_mode not in ('strict','exact_support'):
        raise ValueError('covariance_mode must be strict or exact_support')
    if dps < 30:
        raise ValueError('use at least 30 decimal digits')
    with mp.workdps(dps):
        y = _as_vector(observed)
        directions = [_as_vector(v) for v in basis]
        if not directions or any(v.rows != y.rows for v in directions):
            raise ValueError('need model directions of the observed dimension')
        pinv,kept,condition,values,vectors,roundoff = _spectral(covariance,relative_tolerance)
        if pinv.rows != y.rows:
            raise ValueError('covariance and observation disagree on dimension')
        rank,width = len(kept),len(directions)
        extra = {}
        if rank < y.rows:
            if covariance_mode == 'strict':
                raise ValueError('singular or truncated covariance: exact support must be declared; '
                                 'empirical low rank does not establish deterministic constraints')
            statistic,degrees,amplitudes,cr,mr,st = _exact_support_fit(
                y,directions,kept,values,vectors,mp.mpf(relative_tolerance),roundoff)
            extra = {'covariance_mode':'exact_support','deterministic_constraint_rank':cr,
                     'fitted_stochastic_rank':mr,'support_numerical_tolerance':str(st)}
        else:
            # Preserve the original full-rank arithmetic and return fields.
            gram = mp.zeros(width,width)
            rhs = mp.zeros(width,1)
            weighted = [pinv*v for v in directions]
            for i in range(width):
                for j in range(width):
                    gram[i,j] = (directions[i].T*weighted[j])[0]
                rhs[i] = (directions[i].T*(pinv*y))[0]
            if _rank_of(gram,relative_tolerance) < width:
                raise ValueError('model directions are linearly dependent after weighting; drop one')
            amplitudes = mp.lu_solve(gram,rhs)
            residual = y-sum((amplitudes[i]*directions[i] for i in range(width)),mp.zeros(y.rows,1))
            statistic = (residual.T*(pinv*residual))[0]
            degrees = rank-width
            if degrees < 0:
                raise ValueError('model has more directions than covariance rank')
            if degrees == 0:
                # The design spans all observations. Do not turn roundoff into p=0.
                tolerance = mp.power(10, -mp.dps//2)
                scale = max(1, mp.norm(y), mp.norm(y-residual))
                if mp.norm(residual) > tolerance*scale:
                    raise ArithmeticError('saturated solve lost numerical accuracy')
                statistic = mp.mpf(0)
        return {'statistic':float(statistic),'degrees_of_freedom':degrees,
                'amplitudes':[float(x) for x in amplitudes], 'covariance_rank':rank,
                'covariance_condition_number':float(condition) if condition is not None else None,
                'p_value':chi_square_upper_tail(float(statistic),degrees),
                'equivalent_sigma':_equivalent_sigma(float(statistic),degrees), **extra}


def ray_residual(observed,covariance,direction,
                 relative_tolerance=DEFAULT_RELATIVE_TOLERANCE,dps=DEFAULT_DPS,
                 *,covariance_mode='strict'):
    """Signed-real line test; in two nonsingular dimensions this is Fieller z^2."""
    return subspace_residual(observed,covariance,[direction],relative_tolerance,dps,
                             covariance_mode=covariance_mode)


def _rank_of(matrix,relative_tolerance):
    values,_ = mp.eigsy(matrix)
    magnitudes = [abs(values[i]) for i in range(matrix.rows)]
    largest = max(magnitudes)
    return sum(x > mp.mpf(relative_tolerance)*largest for x in magnitudes) if largest else 0


def chi_square_upper_tail(statistic,degrees):
    if degrees <= 0:
        return 1.0 if statistic <= 0 else 0.0
    if statistic <= 0:
        return 1.0
    with mp.workdps(40):
        return float(mp.gammainc(mp.mpf(degrees)/2,mp.mpf(statistic)/2,mp.inf,regularized=True))


def _equivalent_sigma(statistic,degrees):
    tail = chi_square_upper_tail(statistic,degrees)
    if tail <= 0:
        return float('inf')
    if tail >= 1:
        return 0.0
    with mp.workdps(40):
        return float(-mp.sqrt(2)*mp.erfinv(mp.mpf(tail)-1))
