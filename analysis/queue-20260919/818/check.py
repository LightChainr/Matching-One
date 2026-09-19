"""Exact postprocessing of existing triangular C and integer Walsh energies.
No rank/winding implementation and no new configuration enumeration.
"""
from pathlib import Path
from fractions import Fraction as F
import json
import math

ROOT=Path(__file__).resolve().parent

def main():
    source=json.loads((ROOT/'quoted-input.json').read_text());rows=[]
    for s in source['systems']:
        L=s['L'];n=L*L;c=s['C'];energy=s['walsh_energy'];total=1<<n
        assert len(c)==3 and all(len(a)==n+1 for a in c) and len(energy)==n+1
        assert all(type(v) is int and v>=0 for a in c for v in a)
        assert all(type(v) is int and v>=0 for v in energy)
        assert all(sum(a[k] for a in c)==math.comb(n,k) for k in range(n+1))
        assert all(c[j][k]==c[2-j][n-k] for j in range(3) for k in range(n+1))
        totals=list(map(sum,c));assert totals[0]==totals[2]
        var=F(totals[0]+totals[2],total)
        assert F(sum(energy),total*total)==var
        assert all(energy[k]==0 for k in range(0,n+1,2))
        d=[b-a for a,b in zip(c[0],c[2])]
        mp=F(sum((2*k-n)*v for k,v in enumerate(d)),1<<(n-1))
        ed=mp/n
        ed2=F(4*sum(k*v for k,v in enumerate(energy)),n*total*total)
        assert ed2==F(s['quoted_E_delta_squared'])
        p2=(ed2-ed)/2;p1=2*ed-ed2;p0=1-p1-p2
        probs=[p0,p1,p2]
        assert all(p>=0 for p in probs) and sum(probs)==1
        counts=[p*(1<<(n-1)) for p in probs]
        assert all(a.denominator==1 for a in counts)
        assert p1+2*p2==ed and p1+4*p2==ed2
        rows.append({'L':L,'p':'1/2','N':n,'quoted_configuration_batch':total,
            'rank_counts':totals,'P_rank':[str(F(a,total)) for a in totals],
            'P0_minus_P2_exact':'0','C_columns_checked':n+1,
            'coefficient_complement_identity':True,'walsh_parseval_residual':'0',
            'even_walsh_energy':'0','Mprime':str(mp),'E_delta':str(ed),'E_delta_squared':str(ed2),
            'delta_probabilities':[str(a) for a in probs],
            'derived_fixed_site_delta_counts':[a.numerator for a in counts],
            'R_jump2_exact':str(2*p2/ed),'R_jump2_display':format(float(2*p2/ed),'.12f')})
    out={'issue':818,'status':'exact small-size postprocessing complete',
         'systems':rows,'new_configurations_enumerated':0,'rank_kernel_invoked':False,
         'direct_configuration_duality_rerun':False,
         'source':'Existing PR773 integer rank counts and unnormalized integer Walsh energies',
         'scope':'Integer/Fraction reconstruction, not a rerun of the historical enumeration or a scaling claim.'}
    (ROOT/'result.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(out,ensure_ascii=False,indent=2))

if __name__=='__main__':
    main()
