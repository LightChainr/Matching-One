from __future__ import annotations
import argparse, json
from pathlib import Path
import mpmath as mp


def compute(src):
    rows=[]
    for r in src['finite']:
        L=mp.mpf(r['L']); p=mp.mpf(r['root']); q=1-p; a=mp.mpf(r['a'])
        g1,g2,g3=map(mp.mpf,r['conditional_K_cumulants']['delta_g1_g2_g3'])
        N=L*L
        cov=a*g1
        varx=2*a
        vark=N*p*q
        rho2=cov*cov/(varx*vark)
        z1=-2/g1
        rows.append({
            'L':int(L),
            'g1_over_L_3over4':mp.nstr(g1/L**mp.mpf('.75'),40),
            'g2_over_L_3over4':mp.nstr(g2/L**mp.mpf('.75'),40),
            'g3_over_L_9over4':mp.nstr(g3/L**mp.mpf('2.25'),40),
            'thermal_logit_curvature_candidate_g2_over_g1':mp.nstr(g2/g1,40),
            'metric_free_odd_shape_candidate_g3_over_g1_cubed':mp.nstr(g3/g1**3,40),
            'rho2_K_X':mp.nstr(rho2,40),
            'sqrtL_times_rho2':mp.nstr(mp.sqrt(L)*rho2,40),
            'L_3over4_times_dzds_topological_source':mp.nstr(L**mp.mpf('.75')*z1,40),
            'identities':{
                'varK':'N p(1-p)',
                'varX':'2a',
                'covKX':'a DeltaE[K]',
                'rho2':'p(1-p) (M_p_prime)^2 / (2 a N)',
                'dzds':'-Var(X)/Cov(K,X)=-2/DeltaE[K]'
            }
        })
    return {'scope':'derived exact-finite diagnostics from committed rank-sector histograms; two sizes are controls, not an exponent fit',
            'conditional_hypothesis':'if g_L(z)=G(L^(3/4)t(z))+o(1) with G odd, g1,g2 scale L^(3/4), g3 scales L^(9/4); g2/g1 probes t_zz/t_z and g3/g1^3 is metric-free',
            'rows':rows}

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--input',type=Path,required=True); ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args(); mp.mp.dps=70
    out=compute(json.loads(a.input.read_text())); a.output.write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps(out,indent=2))
