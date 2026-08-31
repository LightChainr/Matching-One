#!/usr/bin/env python3
"""Common-label Euler-invisible thermal curves; reads saved integer hist only."""
import csv
from hashlib import sha256
import json
from pathlib import Path
import subprocess

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import betainc
from scipy.stats import binom

from p334_euler_invisible_thermal_shape import reduced_roots

ROOT = Path(__file__).resolve().parents[1]
SOURCE = "4db356e1b026853468f94d59d938895a2367ceb7"
HISTOGRAM = "results/p334-common-label-euler-tangent/signed_birth_histograms.json"
OUT = ROOT/"results/p334-paired-euler-thermal-response"
CHANNELS = ("plus_S", "plus_D", "minus_S", "minus_D",
            "cross_orientation_antisym", "same_orientation_difference")
OBSERVERS = ("F1", "F2", "A", "E")


def four_fields(two_births):
    return np.stack((two_births[..., 0, :], two_births[..., 1, :],
                     two_births[..., 0, :]+two_births[..., 1, :],
                     two_births[..., 1, :]-two_births[..., 0, :]), axis=-2)


def summary(values):
    return {"mean":values.mean(0).tolist(),
            "se":(values.std(0,ddof=1)/np.sqrt(20)).tolist(),
            "joint20_batch_values":values.tolist()}


def main():
    data = subprocess.check_output(["git","show",f"{SOURCE}:{HISTOGRAM}"],cwd=ROOT)
    source = json.loads(data)
    denominator = source["batch_denominator"]
    pref = source["p_ref"]
    grid = np.unique(np.r_[np.linspace(0,1,1001),pref])
    result = {"histogram_commit":SOURCE,"histogram_path":HISTOGRAM,
              "histogram_sha256":sha256(data).hexdigest(),"batch_denominator":denominator,
              "channels":list(CHANNELS),"observers":list(OBSERVERS),"p_ref":pref,
              "curve_method_source":"7c60b8a7: exact integer cumulative Bernstein coefficients and endpoint-factor root readout",
              "definitions":{
                  "S":"(first+second)/2","D":"(first-second)/original delta_cos4",
                  "A":"H1+H2","E":"H2-H1",
                  "cross_orientation_antisym":"delta_cos4*plus_D/2-minus_S=(T21-T12)/2; no reciprocal-symmetry assumption",
                  "same_orientation_difference":"delta_cos4*plus_D/2+minus_S=(T11-T22)/2",
                  "integral":"Hj integrates to minus the signed Kj first moment divided by N+1",
                  "g_half":"The producer denominator64000 already contains g_plus/minus half factors; no second mark factor"},
              "inference":"All channels and every p share the same20 original paired batches. Pointwise SEs and post-readout extrema are not simultaneous/selection-adjusted intervals. Root labels are descriptive, never a new test or a fitted physical phase.",
              "new_MC":0,"new_histogram_extraction":0,"sizes":{}}
    OUT.mkdir(parents=True,exist_ok=True)
    with (OUT/"thermal_curve.csv").open("w",newline="") as stream:
        writer=csv.writer(stream)
        writer.writerow(["N","channel","p",*[f"mean_{o}" for o in OBSERVERS],*[f"pointwise_se_{o}" for o in OBSERVERS]])
        for ns, src in source["sizes"].items():
            n=int(ns)
            hist=np.asarray(src["batch_integer_histograms"],dtype=np.int64)
            if hist.shape != (20,5,2,2,2,n+1) or np.any(hist.sum(-1)):
                raise ValueError("Expected complete zero-mass signed integer histograms on the fixed axes")
            c=np.cumsum(hist,axis=-1).sum(axis=1)
            # Orientation sum/difference remains integer until the common divisor.
            s=c[:,:,0]+c[:,:,1]
            d=c[:,:,0]-c[:,:,1]
            pair=np.stack((s[:,0],d[:,0],s[:,1],d[:,1],d[:,0]-s[:,1],d[:,0]+s[:,1]),axis=1)
            coef=four_fields(pair)
            delta=src["delta_cos4"]
            divisors=np.array([2,delta,2,delta,2,2])*denominator
            kk=np.arange(n+1)
            basis=binom.pmf(kk[:,None],n,grid[None,:])
            curves=(coef@basis)/divisors[None,:,None,None]
            means=curves.mean(0);ses=curves.std(0,ddof=1)/np.sqrt(20)
            at_pref=(coef@binom.pmf(kk,n,pref))/divisors[None,:,None]
            integrated=coef.sum(-1)/((n+1)*divisors[None,:,None])
            summed_hist=hist.sum(axis=1)
            hs=summed_hist[:,:,0]+summed_hist[:,:,1]
            hd=summed_hist[:,:,0]-summed_hist[:,:,1]
            hpair=np.stack((hs[:,0],hd[:,0],hs[:,1],hd[:,1],hd[:,0]-hs[:,1],hd[:,0]+hs[:,1]),axis=1)
            clock_moments=(hpair@kk)/divisors[None,:,None]
            joint=np.column_stack((at_pref.reshape(20,-1),integrated.reshape(20,-1)))
            centered=joint-joint.mean(0)
            key_labels=[f"{endpoint}.{channel}.{observer}" for endpoint in ("p_ref","integral") for channel in CHANNELS for observer in OBSERVERS]
            size={"batch_ids":src["batch_ids"],"delta_cos4":delta,
                  "key_joint_labels":key_labels,"key_joint20_batch_vectors":joint.tolist(),
                  "key_covariance_of_mean":(centered.T@centered/(20*19)).tolist(),
                  "p_ref":summary(at_pref),"integral":summary(integrated),
                  "signed_K1_K2_moments":summary(clock_moments),"channels":{}}
            for ci, channel in enumerate(CHANNELS):
                for j,p in enumerate(grid):
                    writer.writerow([n,channel,p,*means[ci,:,j],*ses[ci,:,j]])
                item={"p_ref":summary(at_pref[:,ci]),"integral":summary(integrated[:,ci]),"shape":{}}
                for oi in (2,3):
                    coefficients=coef[:,ci,oi]
                    divisor=divisors[ci]
                    aggregate=coefficients.sum(0)
                    roots=reduced_roots(aggregate)

                    def value_at(p):
                        return coefficients@binom.pmf(kk,n,p)/divisor

                    bounds=[0.,*roots["roots"],1.]
                    lobes=[]
                    for left,right in zip(bounds[:-1],bounds[1:]):
                        candidates=np.r_[left,grid[(grid>left)&(grid<right)],right]
                        values=np.array([value_at(p).mean() for p in candidates])
                        peak_index=int(np.argmax(np.abs(values)))
                        sign=1 if values[peak_index]>=0 else -1
                        lower=candidates[max(0,peak_index-1)];upper=candidates[min(len(candidates)-1,peak_index+1)]
                        if upper>lower:
                            opt=minimize_scalar(lambda p:-sign*value_at(p).mean(),bounds=(lower,upper),method="bounded",options={"xatol":2e-12})
                            peak=float(opt.x)
                        else: peak=float(candidates[peak_index])
                        batch_peak=value_at(peak)
                        peak_se=float(batch_peak.std(ddof=1)/np.sqrt(20))
                        peak_mean=float(batch_peak.mean())
                        z=abs(peak_mean)/peak_se if peak_se else None
                        integrated_basis=(betainc(kk+1,n-kk+1,right)-betainc(kk+1,n-kk+1,left))/(n+1)
                        area=coefficients@integrated_basis/divisor
                        lobes.append({"interval":[left,right],"sign":sign,"peak_p":peak,
                                      "peak_mean":peak_mean,"peak_pointwise_se":peak_se,
                                      "peak_abs_z":z,"peak_joint20_values":batch_peak.tolist(),
                                      "area_mean":float(area.mean()),"area_se_at_point_root_limits":float(area.std(ddof=1)/np.sqrt(20)),
                                      "area_joint20_values":area.tolist(),
                                      "birth_response_at_peak":summary((coef[:,ci,:2]@binom.pmf(kk,n,peak))/divisor)})
                    root_details=[]
                    derivative=n*np.diff(coefficients,axis=1)/divisor
                    for j,r in enumerate(roots["roots"]):
                        rvalue=value_at(r)
                        slope=float((derivative@binom.pmf(np.arange(n),n-1,r)).mean())
                        flank_z=[lobes[j]["peak_abs_z"],lobes[j+1]["peak_abs_z"]]
                        resolved=all(z is not None and z>=2 for z in flank_z)
                        root_details.append({"p":r,"flanking_peak_abs_z":flank_z,
                                             "two_flanks_exceed_2_pointwise_SE":resolved,
                                             "slope":slope,
                                             "local_delta_se":float(rvalue.std(ddof=1)/np.sqrt(20)/abs(slope)) if resolved and slope else None,
                                             "interpretation":"Descriptive sign exchange; pointwise amplitudes, not a certified/selection-adjusted root" if resolved else "Weak numerical crossing, not a supported response phase"})
                    item["shape"][OBSERVERS[oi]]={"root_numerics":roots,"root_details":root_details,"lobes":lobes}
                size["channels"][channel]=item
                if channel in ("plus_D","minus_D","plus_S","minus_S"):
                    print(n,channel,"pref A/E",item["p_ref"]["mean"][2:],"SE",item["p_ref"]["se"][2:],
                          "integral A/E",item["integral"]["mean"][2:],flush=True)
                    for obs in ("A","E"):
                        print(obs,[(round(l["peak_p"],6),l["peak_mean"],l["peak_abs_z"]) for l in item["shape"][obs]["lobes"]],flush=True)
            result["sizes"][ns]=size
    result["script_sha256"]=sha256(Path(__file__).read_bytes()).hexdigest()
    (OUT/"score.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    main()
