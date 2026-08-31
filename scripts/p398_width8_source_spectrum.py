#!/usr/bin/env python3
"""P398: width-eight positive join/detach source spectrum, no Monte Carlo."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import scipy
from scipy import linalg, sparse
from scipy.optimize import brentq
from scipy.sparse.linalg import spsolve

from noncrossing_connectivity_codec import canonical_rgs, noncrossing_states
from p321_homology_trace_certificate import join_adjacent, rotate_state
from p333_generic_q_detach_intertwiner import detach_state

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/p398-width8-source-spectrum/latest.json"


def kreweras(state):
    """Noncrossing complement: cycles of p^-1 c, with c=(0 1 ... w-1)."""
    n = len(state)
    permutation = list(range(n))
    for label in set(state):
        block = [j for j in range(n) if state[j] == label]
        for current, following in zip(block, block[1:]+block[:1]):
            permutation[current] = following
    inverse = [permutation.index(j) for j in range(n)]
    complement = [inverse[(j+1) % n] for j in range(n)]
    labels = [-1]*n
    for j in range(n):
        if labels[j] < 0:
            current = j
            while labels[current] < 0:
                labels[current] = j
                current = complement[current]
    return canonical_rgs(labels)


def model(width):
    assert width in (4, 8)
    states = noncrossing_states(width)
    index = {state: j for j, state in enumerate(states)}
    n = len(states)
    rows, cols, values = [], [], []
    for source, state in enumerate(states):
        for site in range(width):
            for action in (join_adjacent, detach_state):
                rows.extend((index[action(state, site)], source))
                cols.extend((source, source))
                values.extend((1, -1))
    forward = sparse.coo_matrix((values, (rows, cols)), shape=(n, n)).tocsc()
    forward.eliminate_zeros()
    assert sparse.csgraph.connected_components(forward,directed=True,connection="strong",return_labels=False) == 1
    mass = -forward.T.tocsr()
    weight = (1j) ** np.arange(width)
    functions = np.array([
        [sum(weight[j] * (state[j] == state[(j+1) % width]) for j in range(width)),
         sum(weight[j] * (state.count(state[j]) == 1) for j in range(width))]
        for state in states
    ])
    size_two = np.array([sum(weight[j] * (state.count(state[j]) == 2)
                             for j in range(width)) for state in states])
    rescue = []
    for state in states:
        value = 0j
        for j in range(width):
            pair = {state[j],state[(j+1) % width]}
            if len(pair) == 2:
                contacts = sum({state[a],state[(a+1) % width]} == pair for a in range(width))
                value += weight[j]*(contacts-1)
        rescue.append(value)
    rescue = np.array(rescue)
    assert np.array_equal(-mass @ functions[:, 1], -3*functions[:, 1]+size_two)
    assert np.array_equal(-mass @ functions[:, 0], -3*functions[:, 0]+rescue)
    # One exact complex degree of freedom for each orbit compatible with character i.
    seen, orbits, representatives = set(), [], []
    for state in states:
        if state in seen:
            continue
        orbit, current = [], state
        while current not in orbit:
            orbit.append(current)
            current = rotate_state(current)
        seen.update(orbit)
        if len(orbit) % 4 == 0:
            orbits.append([index[x] for x in orbit])
            representatives.append(index[state])
    qrows, qcols, qvals, integer_values = [], [], [], []
    for col, orbit in enumerate(orbits):
        for step, row in enumerate(orbit):
            qrows.append(row)
            qcols.append(col)
            qvals.append((1j)**step/np.sqrt(len(orbit)))
            integer_values.append((1j)**step)
    basis = sparse.csc_matrix((qvals, (qrows, qcols)), shape=(n, len(orbits)))
    integer_basis = sparse.csc_matrix((integer_values, (qrows, qcols)), shape=basis.shape)
    reduced = (basis.conj().T @ mass @ basis).toarray()
    source = np.asarray(basis.conj().T @ functions)
    assert np.max(np.abs(basis @ source-functions)) < 1e-12
    assert np.max(np.abs(mass @ basis-basis @ reduced)) < 1e-12
    reduced_integer = (mass @ integer_basis)[representatives, :].toarray()
    source_integer = functions[representatives]
    complement = np.array([index[kreweras(state)] for state in states])
    assert (mass[complement, :][:, complement]-mass).nnz == 0
    assert np.array_equal(functions[complement, 0], -1j*functions[:, 1])
    assert np.array_equal(functions[complement, 1], functions[:, 0])
    assert all(states[complement[complement[j]]] == rotate_state(state, -1)
               for j, state in enumerate(states))
    complement_reduced = integer_basis[complement[representatives], :].toarray()
    assert np.trace(complement_reduced) == 0
    equations = forward.tolil().astype(float)
    equations[-1, :] = np.ones(n)
    rhs = np.zeros(n)
    rhs[-1] = 1
    stationary = spsolve(equations.tocsc(), rhs)
    assert np.min(stationary) > 0 and np.max(np.abs(forward @ stationary)) < 1e-12
    assert np.max(np.abs(stationary @ functions)) < 1e-12
    return states, mass, functions, size_two, basis, reduced, source, reduced_integer, source_integer, stationary


def modular_krylov(matrix, source, ray_sign=None):
    """Rank lower bound over Q(i), by exact reduction i=256 mod prime 65537."""
    prime, root = 65537, 256
    convert = lambda value: (np.rint(value.real).astype(np.int64)
                             + root*np.rint(value.imag).astype(np.int64)) % prime
    operator, block = convert(matrix), convert(source)
    if ray_sign is not None:
        sqrt_two = next(x for x in range(prime) if x*x % prime == 2)
        phase = (1-root)*pow(sqrt_two, -1, prime) % prime
        assert phase*phase % prime == -root % prime
        block = (block[:, :1]+ray_sign*phase*block[:, 1:2]) % prime
    basis, pivots, history = [], [], []
    for power in range(matrix.shape[0]+1):
        for column in block.T:
            vector = column.copy()
            for pivot, previous in zip(pivots, basis):
                vector = (vector-vector[pivot]*previous) % prime
            nonzero = np.flatnonzero(vector)
            if len(nonzero):
                pivot = int(nonzero[0])
                vector = vector*pow(int(vector[pivot]), -1, prime) % prime
                pivots.append(pivot)
                basis.append(vector)
        history.append(len(basis))
        if len(basis) == matrix.shape[0] or (power and history[-1] == history[-2]):
            break
        block = operator @ block % prime
    return {"prime":prime,"i_image":root,"rank_by_max_power":history,
            "reachable_rank_lower_bound":len(basis),
            "ray_sign":ray_sign,
            "full_sector_reached":len(basis)==matrix.shape[0]}


def complex_display(value):
    value = np.asarray(value)
    return np.stack((value.real, value.imag), axis=-1).tolist()


def fingerprint(matrix):
    trace, determinant = np.trace(matrix).real, np.linalg.det(matrix).real
    split = 2*np.arccosh(trace/(2*np.sqrt(determinant)))/(-np.log(determinant))
    mixing = (4*matrix[0,1]*matrix[1,0]/(trace*trace-4*determinant)).real
    return {"I_mass":float(split),"I_mix":float(mixing),
            "linked_invariant":float(25*split*split*(2-mixing))}


def analyze(width):
    states, mass, f, t2, q, h, source, hi, fi, pi = model(width)
    covariance = f.T @ (pi[:,None]*f.conj())
    propagated = mass @ f
    best = np.linalg.solve(f.conj().T @ (pi[:,None]*f), f.conj().T @ (pi[:,None]*propagated))
    leakage = propagated-f @ best
    leakage_covariance = leakage.conj().T @ (pi[:,None]*leakage)
    eigenvalues, eigenvectors = linalg.eig(h)
    order = np.lexsort((eigenvalues.imag, eigenvalues.real))
    eigenvalues, eigenvectors = eigenvalues[order], eigenvectors[:,order]
    source_coefficients = linalg.solve(eigenvectors, source)
    full_modes = q @ eigenvectors
    # C(s)=sum_j residue_j exp(-conjugate(m_j)*s), from the same physical pi.
    left = f.T @ (pi[:,None]*full_modes.conj())
    residues = np.array([np.outer(left[:,j],source_coefficients[j,:].conj())
                         for j in range(len(eigenvalues))])
    assert np.max(np.abs(sum(residues)-covariance)) < 1e-9
    ray_change = np.array([[1,1],[-np.exp(-1j*np.pi/4),np.exp(-1j*np.pi/4)]])/np.sqrt(2)
    ray_residues = np.array([ray_change.T @ row @ ray_change.conj() for row in residues])
    ray_c0 = ray_change.T @ covariance @ ray_change.conj()
    assert np.max(np.abs(ray_residues[:,0,1])) < 1e-9
    assert np.max(np.abs(ray_residues[:,1,0])) < 1e-9
    ray_masses = [[],[]]
    for j, value in enumerate(eigenvalues):
        ray = int(np.argmax(np.abs(np.diag(ray_residues[j]))))
        ray_masses[ray].append({"mass":float(value.real),
                               "spectral_residue":float(ray_residues[j,ray,ray].real)})
    low_modes=[]
    for j in range(min(16,len(eigenvalues))):
        low_modes.append({"mass_re_im":complex_display(eigenvalues[j]),
            "residue_re_im":complex_display(residues[j]),
            "psi_minus_plus_residue_re_im":complex_display(ray_residues[j]),
            "residue_frobenius_norm":float(linalg.norm(residues[j]))})
    samples=[]
    kernels={}
    for distance in (.05,.1,.25,.5,1.,2.):
        corr=np.sum(residues*np.exp(-eigenvalues.conj()*distance)[:,None,None],axis=0)
        u=linalg.solve(covariance,corr)
        kernels[distance]=u
        ray_corr=ray_change.T @ corr @ ray_change.conj()
        samples.append({"s":distance,"C_re_im":complex_display(corr),
                        "U_re_im":complex_display(u),
                        "psi_minus_plus_covariance":np.diag(ray_corr).real.tolist(),
                        **fingerprint(u)})
    def normalized_ray_difference(distance):
        decay=np.exp(-eigenvalues.conj()*distance)
        ray_corr=np.sum(ray_residues*decay[:,None,None],axis=0)
        return float((ray_corr[0,0]/ray_c0[0,0]-ray_corr[1,1]/ray_c0[1,1]).real)
    crossing=None
    if width==8:
        crossing=float(brentq(normalized_ray_difference,.1,.5))
    b4=np.array([[5,1+1j],[1-1j,5]])
    defect=propagated-f@b4
    pair_groups={}
    witness=None
    for j in range(len(states)):
        key=tuple(f[j])
        if key in pair_groups:
            other=pair_groups[key]
            if not np.array_equal(propagated[j],propagated[other]):
                witness={"states":[states[other],states[j]],"same_A_L":complex_display(f[j]),
                         "mass_A_L":complex_display(propagated[[other,j]]),
                         "size_two_charge":complex_display(t2[[other,j]])}
                break
        else:
            pair_groups[key]=j
    return {"width":width,"states":len(states),"cyclic_character":"i",
        "lattice_momentum_k":width//4,"exact_character_sector_dimension":h.shape[0],
        "source_krylov":modular_krylov(hi,fi),
        "kreweras":{"permutation_commutes_with_G_exact":True,
            "squared_is_rotation_minus_one_exact":True,
            "readout_action_exact":"A(K state)=-i L(state); L(K state)=A(state)",
            "trace_on_character_i_exact":0,
            "exact_minus_plus_dimensions":[h.shape[0]//2,h.shape[0]//2],
            "source_minus_krylov":modular_krylov(hi,fi,-1),
            "source_plus_krylov":modular_krylov(hi,fi,1)},
        "stationary_min":float(pi.min()),"stationary_max":float(pi.max()),
        "C0_re_im":complex_display(covariance),
        "best_stationary_L2_two_channel_mass_re_im":complex_display(best),
        "two_channel_leakage_covariance_re_im":complex_display(leakage_covariance),
        "first_derivative_new_directions":"R=weighted extra adjacent contacts between the two endpoint clusters; T2=weighted size-two cluster membership",
        "width4_identity_max_defect":float(np.max(np.abs(defect))),
        "same_readouts_different_next_derivative_witness":witness,
        "spectrum_max_imaginary_part":float(np.max(np.abs(eigenvalues.imag))),
        "all_masses_re_im":complex_display(eigenvalues),"lowest_modes":low_modes,
        "psi_minus_plus_C0":np.diag(ray_c0).real.tolist(),
        "psi_minus_plus_spectrum":ray_masses,
        "two_lowest_mass_ratio":float(eigenvalues[1].real/eigenvalues[0].real),
        "nonzero_normalized_ray_crossing_s":crossing,
        "minimum_spectral_residue":float(np.min(np.diagonal(ray_residues,axis1=1,axis2=2).real)),
        "max_eigenpair_residual":float(np.max(np.abs(h@eigenvectors-eigenvectors*eigenvalues))),
        "max_cross_ray_spectral_residue":float(max(np.max(abs(ray_residues[:,0,1])),np.max(abs(ray_residues[:,1,0])))),
        "kernel_samples":samples,
        "semigroup_relative_defect_at_half_plus_half":float(linalg.norm(kernels[1.]-kernels[.5]@kernels[.5])/linalg.norm(kernels[1.])),
        "eigenvector_condition_number":float(np.linalg.cond(eigenvectors))}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json",type=Path,default=OUT)
    args=parser.parse_args()
    value={"schema":"matching-one/p398-width8-source-spectrum/v1",
        "parent":"dbd408154b4215ca41fbf26c0fd962997074f05d",
        "generator":"G=sum_j(J_j-I)+sum_j(D_j-I), Q=1, continuous self-dual limit",
        "readouts":{"A":"sum_j i^j 1{j connected to j+1}",
                    "L":"sum_j i^j 1{j is singleton}",
                    "T2":"sum_j i^j 1{block of j has size 2}",
                    "R":"sum_j i^j 1{j,j+1 disconnected}*(number of adjacent edges connecting their two blocks minus 1)",
                    "exact_identity":"G L=-3L+T2; G A=-3A+R"},
        "arithmetic":"integer sparse generator, exact orbit dimension and modular rank lower bounds; float64 stationary/spectral display",
        "environment":{"numpy":np.__version__,"scipy":scipy.__version__},
        "rows":[]}
    for width in (4,8):
        row=analyze(width)
        value["rows"].append(row)
        print(width,row["exact_character_sector_dimension"],row["source_krylov"],flush=True)
        print("low masses",[x["mass_re_im"] for x in row["lowest_modes"][:8]],flush=True)
        print("fingerprints",[(x["s"],x["linked_invariant"]) for x in row["kernel_samples"]],flush=True)
        print("leakage",row["two_channel_leakage_covariance_re_im"],flush=True)
    paths=[Path(__file__),*(ROOT/"scripts"/name for name in (
        "noncrossing_connectivity_codec.py","p321_homology_trace_certificate.py","p333_generic_q_detach_intertwiner.py"))]
    value["input_sha256"]={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    args.json.parent.mkdir(parents=True,exist_ok=True)
    args.json.write_text(json.dumps(value,indent=2)+"\n")


if __name__=="__main__":
    main()
