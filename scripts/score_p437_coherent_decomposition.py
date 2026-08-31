#!/usr/bin/env python3
"""Archive-only fixed-S degree-five coefficient versus higher-degree energy."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

from integer_period_torus import integer_torus_geometry, matrix_vector
from p437_fixed_support_mc import SUPPORT, DENOMINATOR, energy_numerator
from p437_high_pass_mc import MATRICES
from p437_positive_difference_bridge import H5, spectral_multiplier

D4 = (((1,0),(0,1)), ((0,-1),(1,0)), ((-1,0),(0,-1)), ((0,1),(-1,0)),
      ((-1,0),(0,1)), ((1,0),(0,-1)), ((0,1),(1,0)), ((0,-1),(-1,0)))


def maps_lattice(source, target, rotation):
    matrix = source.periods.matrix
    try:
        for column in ((matrix[0][0],matrix[1][0]), (matrix[0][1],matrix[1][1])):
            target.periods.winding(matrix_vector(rotation,column))
    except ValueError:
        return False
    return True  # equal index 112: lattice containment is equality


def transported_support(source, target, rotation, translation, dual=False):
    output = []
    for index in SUPPORT:
        edge = source.primal_edges[index]
        origin = source.coordinates[edge.i]
        direction = (edge.dx,edge.dy)
        if dual:
            if direction == (1,0):
                origin, direction = (origin[0],origin[1]-1), (0,1)
            else:
                origin, direction = (origin[0]-1,origin[1]), (1,0)
        origin, direction = matrix_vector(rotation,origin), matrix_vector(rotation,direction)
        origin = (origin[0]+translation[0],origin[1]+translation[1])
        if direction in ((-1,0),(0,-1)):
            origin = (origin[0]+direction[0],origin[1]+direction[1])
            direction = (-direction[0],-direction[1])
        output.append(2*target.vertex(origin)+int(direction==(0,1)))
    return sorted(output)


def symmetry_certificate():
    geometries = [integer_torus_geometry(matrix) for matrix in MATRICES]
    relations = []
    for dual in (False,True):
        for i, source in enumerate(geometries):
            for j, target in enumerate(geometries):
                witnesses = []
                for rotation in D4:
                    if not maps_lattice(source,target,rotation):
                        continue
                    for translation in target.coordinates:
                        if transported_support(source,target,rotation,translation,dual)==list(SUPPORT):
                            witnesses.append({"rotation":rotation,"translation":translation})
                if witnesses:
                    relations.append({"source_child":i,"target_child":j,"coefficient_sign":-1 if dual else 1,
                                      "dual_complement":dual,"number_of_maps":len(witnesses),"first_witness":witnesses[0]})
    nontrivial = [r for r in relations if r["source_child"]!=r["target_child"] or r["coefficient_sign"]==-1]
    if [(r["source_child"],r["target_child"],r["coefficient_sign"]) for r in nontrivial] != [(1,2,1),(2,1,1)]:
        raise AssertionError("fixed-support symmetry relation changed")
    return {"scope":"all translations and physical D4 maps between these period lattices, with/without geometric dual complement",
            "support":SUPPORT,"relations":relations,
            "independent_constraint_on_real_child_coefficients":[[0,1,-1]],
            "allowed_child_coefficient_basis":[[1,0,0],[0,1,1]],
            "conclusion":"a1=a2 exactly, hence Im Fhat(S)=0 and Re Fhat(S)=(a0-a1)/3 is allowed; no anti-invariant fixed-support map was found",
            "boundary":"marginal Fourier-coefficient relation; not a pointwise common-stream conjugation identity or proof from a nonsignificant p value"}


def coherent_weight(z):
    m = len(z)
    total = z.sum(axis=0)
    return float((total@total - np.sum(z*z))/(m*(m-1)))


def parameters(x):
    B = float(x[:,2].mean())
    coherent = coherent_weight(x[:,:2])
    return np.array([B,coherent,B-coherent,coherent/B])


def score(directory):
    metadata = json.loads((directory/"run.json").read_text())
    path = directory/"batches.json"
    if hashlib.sha256(path.read_bytes()).hexdigest()!=metadata["batch_sha256"]:
        raise ValueError("input hash mismatch")
    rows = json.loads(path.read_text())
    vectors = []
    for row in rows:
        signed = np.zeros(3)
        energy = 0
        for value in row["classes"]:
            n = np.array(value["child_difference_numerators"])
            signed += n*value["count"]
            energy += energy_numerator(n)*value["count"]
        d = signed/(32*row["samples"])
        vectors.append([(2*d[0]-d[1]-d[2])/6,np.sqrt(3)*(d[2]-d[1])/6,
                        energy/(DENOMINATOR*row["samples"])])
    x = np.array(vectors)
    m = len(x)
    mean = x[:,:2].mean(axis=0)
    covariance = np.cov(x[:,:2],rowvar=False,ddof=1)/m
    estimate = parameters(x)
    plug_in = float(mean@mean)
    correction = float(np.trace(covariance))
    if not np.isclose(estimate[1],plug_in-correction,rtol=1e-12,atol=1e-25):
        raise AssertionError("cross-batch U-statistic differs from covariance-subtracted norm")
    leave_one = np.array([parameters(np.delete(x,b,axis=0)) for b in range(m)])
    jackknife_cov = np.cov(leave_one,rowvar=False,ddof=1)*(m-1)**2/m
    se = np.sqrt(np.diag(jackknife_cov))
    phase_gradient = np.array([-mean[1],mean[0]])/plug_in
    h5,h6 = float(H5),float(spectral_multiplier(6))
    bound_gradient = np.array([h6,h5-h6])
    refined = float(bound_gradient@estimate[:2])
    refined_se = float(np.sqrt(bound_gradient@jackknife_cov[:2,:2]@bound_gradient))
    return {"schema":"matching-one/p437-fixed-S-coherent-decomposition/v1",
            "source_result_commit":"386db0a","source_batch_sha256":metadata["batch_sha256"],
            "new_samples":0,"batches":m,"samples":sum(r["samples"] for r in rows),"support":SUPPORT,
            "fourier_coefficient":{
                "identity":"E[D_SF]=Fhat(S), the exact degree-five coefficient on this one fixed support",
                "mean_re_im":mean.tolist(),"covariance_of_mean_2x2":covariance.tolist(),
                "se_re_im":np.sqrt(np.diag(covariance)).tolist(),
                "z_re_im":(mean/np.sqrt(np.diag(covariance))).tolist(),
                "phase_degrees":float(np.degrees(np.arctan2(mean[1],mean[0]))),
                "phase_delta_se_degrees":float(np.degrees(np.sqrt(phase_gradient@covariance@phase_gradient)))},
            "energy_decomposition":{
                "parameter_order":["B_S","coherent_degree5_weight","outside_dependent_degree6plus_weight","coherent_fraction"],
                "estimates":estimate.tolist(),"jackknife_se":se.tolist(),"joint_jackknife_covariance":jackknife_cov.tolist(),
                "plug_in_squared_mean":plug_in,"mean_variance_bias_removed":correction,
                "coherent_estimator":"sum_{b!=c} Re(conj(Z_b) Z_c)/(m(m-1)) = |mean Z|²-tr(Cov(mean Z))",
                "higher_order_identity":"B_S-|Fhat(S)|²=sum_{T strictly contains S}|Fhat(T)|²",
                "negative_weight_policy":"keep negative unbiased estimates and label unresolved, never interpret them as negative population energy",
                "fraction_boundary":"ratio is descriptive with delete-one-batch uncertainty, not an exactly unbiased ratio"},
            "symmetry":symmetry_certificate(),
            "secondary_refined_population_lower_bound_parameter":{
                "exact_inequality":"A_HP>=h5 |mu|²+h6 (B_S-|mu|²)",
                "h5":str(H5),"h6":str(spectral_multiplier(6)),"h6_over_h5":"63/32",
                "reason":"The outside-dependent residual contains only Fourier supports strictly larger than S, hence degree>=6",
                "estimate":refined,"batch_jackknife_se":refined_se,
                "old_h5_B_S_estimate":h5*estimate[0],
                "improvement_ratio_estimate":refined/(h5*estimate[0]),
                "improvement_ratio_delta_se":(h6/h5-1)*se[3],
                "inference":"secondary reuse of the original fixed-support block; estimated RHS with uncertainty, not a statistically certain numerical lower bound"},
            "interpretation":"The fifteen-SE localized energy mostly measures outside-dependent higher-degree fluctuations; it must not be relabeled pure fifth-degree energy.",
            "dependency_group":"same p437-N112-fixed-S5-lower-bound-fresh20k-20260831, not new evidence block"}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory",type=Path)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(score(args.directory),indent=2,sort_keys=True,allow_nan=False)+"\n")


if __name__=="__main__":
    main()
