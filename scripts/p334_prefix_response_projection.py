#!/usr/bin/env python3
"""Named fixed-rank-cell feature projections of conditional birth responses.

Predictor clocks and response means are latent prefix quantities. Different
quartets estimate their products without the shared-tail self-product bias.
Exact census features have no conditional Monte Carlo measurement error.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE = "375cd3a12b2b7a87d79148a59f62b95898f9e471"
ORIS = ("first", "second")
CELLS = (("00", "01", "02"), ("00", "10", "20"))
PREDICTORS = ("joint_safe_mass", "own_score_energy", "own_safe_degree", "own_safe_loop", "mu_C", "mu_W")
RESPONSES = ("source_first.C", "source_first.W", "source_second.C", "source_second.W")
VARIABLES = PREDICTORS+RESPONSES
MODELS = {"strength": (1,), "contact": (0, 1, 2, 3), "clock": (4, 5), "contact_clock": tuple(range(6))}
TRI = np.triu_indices(len(VARIABLES))


def quartet_variables(b, h, features):
    """features: (prefix,orientation,4); output: (prefix,quartet,ori,10)."""
    clock = np.stack(((b[..., 0]+b[..., 1])/2, b[..., 1]-b[..., 0]), axis=-1)
    # h axes: prefix,quartet,mark,receiver,feature. H_Lf=H_plus+H_minus.
    physical = np.stack((h[:, :, 0]+h[:, :, 1], h[:, :, 0]-h[:, :, 1]), axis=2)
    response = np.stack(((physical[..., 0]+physical[..., 1])/2,
                          physical[..., 1]-physical[..., 0]), axis=-1)
    response = response.transpose(0, 1, 3, 2, 4).reshape(len(b), b.shape[1], 2, 4)
    f = np.broadcast_to(features[:, None, :, :], (*clock.shape[:-1], 4))
    return np.concatenate((f, clock, response), axis=-1)


def sufficient_batches(z, batch, rankcell):
    q = z.shape[1]
    sm = z.sum(axis=1)
    mean = sm/q
    # Ordered distinct-quartet product, symmetrized automatically by the sum.
    u = (np.einsum("noi,noj->noij", sm, sm)-np.einsum("nqoi,nqoj->noij", z, z))/(q*(q-1))
    diag = np.einsum("noi,noj->noij", mean, mean)
    fields = {}
    for oi, ori in enumerate(ORIS):
        for cell in CELLS[oi]:
            mask = rankcell == 3*int(cell[0])+int(cell[1])
            fields[f"{ori}.{cell}.mass"] = mask.astype(float)
            for j, v in enumerate(VARIABLES):
                fields[f"{ori}.{cell}.mean.{v}"] = mask*mean[:, oi, j]
            for i, j in zip(*TRI):
                name = f"{VARIABLES[i]}|{VARIABLES[j]}"
                fields[f"{ori}.{cell}.U.{name}"] = mask*u[:, oi, i, j]
                fields[f"{ori}.{cell}.diag.{name}"] = mask*diag[:, oi, i, j]
    x = np.column_stack(list(fields.values()))
    ids = np.unique(batch)
    counts = np.array([(batch == v).sum() for v in ids])
    if not np.all(counts == counts[0]):
        raise ValueError("Original prefix batches must have equal sizes")
    return list(fields), np.array([x[batch == v].mean(axis=0) for v in ids]), ids, counts


def from_triangle(d, stem):
    out = np.zeros((len(VARIABLES), len(VARIABLES)))
    for i, j in zip(*TRI):
        out[i, j] = out[j, i] = d[stem+f"{VARIABLES[i]}|{VARIABLES[j]}"]
    return out


def solve_projection(k, v):
    scale = np.sqrt(np.diag(k))
    if not np.all(np.isfinite(scale)) or np.any(scale <= 0):
        raise ValueError("A named predictor variance is not identified; no PSD repair applied")
    cor = k/scale[:, None]/scale[None, :]
    eig = np.linalg.eigvalsh(cor)
    if eig[0] <= 0:
        raise ValueError("A named latent predictor Gram is not positive definite; no ridge applied")
    beta = np.linalg.solve(cor, v/scale[:, None])/scale[:, None]
    return beta, eig[0], eig[-1]/eig[0]


def derive(row, labels, population):
    d, out = dict(zip(labels, row)), {}
    for oi, ori in enumerate(ORIS):
        cov = np.zeros((len(VARIABLES), len(VARIABLES)))
        mass = 0.
        for cell in CELLS[oi]:
            stem = f"{ori}.{cell}."
            pi = d[stem+"mass"]; m = int(round(population*pi))
            mu = np.array([d[stem+"mean."+v]/pi for v in VARIABLES])
            u, diag = from_triangle(d, stem+"U."), from_triangle(d, stem+"diag.")
            product = (m*pi*np.outer(mu, mu)-diag)/(m-1)
            cov += u-product
            mass += pi
        out[f"{ori}.receiver_R0_mass"] = mass
        for i, j in zip(*TRI):
            out[f"{ori}.within_G_cov.{VARIABLES[i]}|{VARIABLES[j]}"] = cov[i, j]
        k, v, rr = cov[:6, :6], cov[:6, 6:], cov[6:, 6:]
        own = (2*oi, 2*oi+1)
        observed_loading = 2*v[4, own[0]]-.5*v[5, own[1]]
        out[f"{ori}.own_source.intrinsic_cov_loading"] = observed_loading
        projections = {}
        for model, cols in MODELS.items():
            ix = np.array(cols)
            beta, eig, condition = solve_projection(k[np.ix_(ix, ix)], v[ix])
            proj_var = v[ix].T@beta
            clock_loading = k[4:6, ix]@beta
            projections[model] = proj_var
            out[f"{ori}.{model}.scaled_Gram_min_eigenvalue"] = eig
            out[f"{ori}.{model}.scaled_Gram_condition"] = condition
            for j, response in enumerate(RESPONSES):
                out[f"{ori}.{model}.{response}.projected_response_variance"] = proj_var[j, j]
                out[f"{ori}.{model}.{response}.response_variance_moment"] = rr[j, j]
                for i, col in enumerate(cols):
                    out[f"{ori}.{model}.{response}.beta.{PREDICTORS[col]}"] = beta[i, j]
            if model in ("strength", "contact"):
                predicted = 2*clock_loading[0, own[0]]-.5*clock_loading[1, own[1]]
                out[f"{ori}.{model}.own_intrinsic_cov_loading"] = predicted
                out[f"{ori}.{model}.own_intrinsic_cov_loading_residual"] = observed_loading-predicted
                out[f"{ori}.{model}.own_intrinsic_cov_loading_share"] = predicted/observed_loading
        # Additional geometry information after the two baseline mean clocks.
        for j, response in enumerate(RESPONSES):
            out[f"{ori}.{response}.contact_increment_after_clock"] = projections["contact_clock"][j,j]-projections["clock"][j,j]
            out[f"{ori}.{response}.clock_increment_after_contact"] = projections["contact_clock"][j,j]-projections["contact"][j,j]
        bclock, _, _ = solve_projection(k[4:6, 4:6], v[4:6])
        residual = v[:4]-k[:4, 4:6]@bclock
        for i, feature in enumerate(PREDICTORS[:4]):
            for j, response in enumerate(RESPONSES):
                out[f"{ori}.clock_partial_cov.{feature}|{response}"] = residual[i,j]
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--moments", type=Path, required=True)
    parser.add_argument("--descriptors", type=Path, required=True)
    parser.add_argument("--descriptor-commit", required=True)
    parser.add_argument("--output", type=Path, default=ROOT/"results/p334-prefix-response-projection")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    result = {"schema": "matching-one/p334-prefix-response-projection/v1", "moment_source_commit": SOURCE,
              "descriptor_commit": args.descriptor_commit, "predictors": PREDICTORS, "responses": RESPONSES,
              "models": MODELS, "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "reader_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
              "estimand": "Pooled within-rank-cell latent prefix covariances, zero padded to the full population, restricted to receiver rank0. Separate physical receiver and source, with rank-cell-specific intercepts and common slopes.",
              "boundary": "Exploratory linear projection, not causal attribution or out-of-sample prediction. Clock-only projection reproduces clock-response covariance by definition and is not a closure discovery. Projected variance estimates are not R-squared when latent response variance is noisy.",
              "new_samples": 0, "new_DP": 0, "input_sha256": {}, "sizes": {}}
    for n in (325, 425):
        moment_path, descriptor_path = args.moments/f"N{n}.npz", args.descriptors/f"N{n}.npz"
        with np.load(moment_path, allow_pickle=False) as z, np.load(descriptor_path, allow_pickle=False) as f:
            if not np.array_equal(z["counter"], f["counter"]) or not np.array_equal(z["batch"], f["batch"]):
                raise ValueError("Descriptor and fork prefix identities differ")
            # The descriptor extractor provides these four named columns in
            # physical orientation order; no selection on response values.
            features = np.stack((f["features"][:, [0, 6, 2, 4]],
                                 f["features"][:, [0, 8, 3, 5]]), axis=1)
            a = quartet_variables(z["b"], z["h"], features)
            labels, raw, ids, counts = sufficient_batches(a, z["batch"], z["rankcell"])
        m = len(ids); p = int(counts.sum()); mean = raw.mean(axis=0)
        point = derive(mean, labels, p)
        loo = np.array([list(derive((m*mean-row)/(m-1), labels, p-int(counts[j])).values())
                        for j, row in enumerate(raw)])
        factor = np.sqrt((m-1)/m)*(loo-loo.mean(axis=0))
        result["sizes"][str(n)] = {"batch_ids": ids.tolist(), "prefix_counts": counts.tolist(),
            "raw_labels": labels, "raw_batch_means": raw.tolist(), "labels": list(point),
            "estimate": list(point.values()), "se": np.linalg.norm(factor, axis=0).tolist(),
            "LOO": loo.tolist(), "factor": factor.tolist()}
        for ori in ORIS:
            for suffix in ("own_source.intrinsic_cov_loading", "strength.own_intrinsic_cov_loading_share",
                           "contact.own_intrinsic_cov_loading_share", "contact.own_intrinsic_cov_loading_residual"):
                key = ori+"."+suffix; j = list(point).index(key)
                print(n, key, f"{point[key]:.10g} +/- {np.linalg.norm(factor[:,j]):.6g}", flush=True)
        for path in (moment_path, descriptor_path):
            result["input_sha256"][str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    (args.output/"score.json").write_text(json.dumps(result, separators=(",", ":"), allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
