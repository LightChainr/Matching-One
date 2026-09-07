#!/usr/bin/env python3
"""Exploratory normalized-shape pilot after #610's chart correction.

One attractive scalar fit must not stand in for a full-law fit. Compare both
with identical deciles and covariance inputs. These are conditional plug-in
diagnostics, not calibrated p values or a newly identified critical exponent.
"""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar

ROOT=Path(__file__).resolve().parents[1]


def fit_power(ns, ys, covs):
    dimension=len(ys[0]); whiteners=[]
    for cov in covs:
        value,vector=np.linalg.eigh(cov)
        assert min(value)>0
        whiteners.append((vector/np.sqrt(value)).T)
    def profile(delta):
        design=np.vstack([w@np.column_stack([np.eye(dimension),n**(-delta)*np.eye(dimension)])
                          for n,w in zip(ns,whiteners)])
        obs=np.concatenate([w@y for w,y in zip(whiteners,ys)])
        coefficients=np.linalg.lstsq(design,obs,rcond=None)[0]
        residual=obs-design@coefficients
        return float(residual@residual),coefficients
    opt=minimize_scalar(lambda d:profile(d)[0],bounds=(.05,2),method='bounded')
    chi,coefficients=profile(opt.x)
    return {'delta':float(opt.x),'conditional_chi_square':chi,
            'nominal_residual_dimension':len(ns)*dimension-2*dimension-1,
            'intercept':coefficients[:dimension].tolist(),
            'power_coefficient':coefficients[dimension:].tolist()}


def main():
    data=json.loads((ROOT/'results/astra610-independent/projection-check.json').read_text())['reconstruction']
    g=np.array(json.loads((ROOT/'results/p582-amplitude-law/latest.json').read_text())['frozen_direction'])
    ns=np.array(sorted(map(int,data)),dtype=float);out={}
    for weighting in ('spin0','equal'):
        ref=np.array(data['145'][weighting]['quantiles']);width=ref[7]-ref[1]
        direction=(g-g[4])/width-(ref-ref[4])*(g[7]-g[1])/width**2
        direction/=np.linalg.norm(direction)
        full_y=[];full_s=[];scalar_y=[];scalar_s=[];widths=[]
        indices=[0,1,2,3,5,6,8]  # Q50=0 and Q80-Q20=1 remove two coordinates.
        for n in ns:
            row=data[str(int(n))][weighting];q=np.array(row['quantiles']);w=q[7]-q[1]
            z=(q-q[4])/w
            jac=(np.eye(9)-np.outer(np.ones(9),np.eye(9)[4]))/w-np.outer(q-q[4],np.eye(9)[7]-np.eye(9)[1])/w**2
            cov=jac@row['covariance']@jac.T
            full_y.append(z[indices]);full_s.append(cov[np.ix_(indices,indices)])
            scalar_y.append(np.array([direction@z]));scalar_s.append(np.array([[direction@cov@direction]]))
            widths.append(w)
        out[weighting]={'scalar_direction':direction.tolist(),
                        'scalar_fit':fit_power(ns,scalar_y,scalar_s),
                        'full_seven_coordinate_fit':fit_power(ns,full_y,full_s),
                        'unweighted_width_exponent':float(-np.polyfit(np.log(ns),np.log(widths),1)[0])}
    result={'status':'C2 exploratory; delta-method within-size covariance, cross-size terms omitted; direction learned from same data; no inferential p values',
            'normalization':'Z_N(u)=(Q_N(u)-Q_N(.5))/(Q_N(.8)-Q_N(.2))',
            'model':'Z_N=z_infinity+N^(-delta)*b; full vector fits 7 intercepts+7 loadings+delta',
            'sizes':ns.astype(int).tolist(),'weightings':out}
    (ROOT/'results/astra610-independent/intrinsic-shape-pilot.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
