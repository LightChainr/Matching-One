"""No enumeration and no fitting. Compare published atlas aggregates to rank marginals."""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal, localcontext
import json
import math

ROOT=Path(__file__).resolve().parent

def display(x, digits=22):
    with localcontext() as ctx:
        ctx.prec=digits
        return str(Decimal(x.numerator)/Decimal(x.denominator))

def falling(n,k):
    return math.prod(range(n-k+1,n+1)) if k else 1

def derivative(c,p,order):
    n=len(c)-1;q=1-p;ans=F(0)
    for k,coeff in enumerate(c):
        for a in range(order+1):
            b=order-a
            if a<=k and b<=n-k:
                ans+=coeff*math.comb(order,a)*falling(k,a)*falling(n-k,b)*(-1)**b*p**(k-a)*q**(n-k-b)
    return ans

def main():
    source=json.loads((ROOT/'quoted-input.json').read_text());p=F(source['p']);results=[];ratios=[]
    for s in source['systems']:
        L=s['L'];n=L*L;c=s['C'];d=[b-a for a,b in zip(c[0],c[2])]
        assert all(sum(row[k] for row in c)==math.comb(n,k) for k in range(n+1))
        mp=derivative(d,p,1);mpp=derivative(d,p,2)
        E=F(s['atlas_E_delta']);P2=F(s['atlas_P_jump2']);NJ=F(s['atlas_N_sum_J'])
        r=2*P2/E;res1=n*E-mp;res2=NJ-mpp
        assert res1==res2==0,(L,res1,res2)
        assert 0<=P2<=E/2 and 0<=r<=1
        assert abs(r-F(s['quoted_R_display']))<F(1,10**19)
        ratios.append((L,r))
        results.append({'L':L,'R_jump2':display(r),'R_jump2_exact':str(r),
            'P_jump2':display(P2),'E_delta':display(E),'Mprime_from_C':display(mp),
            'Mdoubleprime_from_C':display(mpp),'N_E_minus_Mprime_exact':str(res1),
            'N_sumJ_minus_Mdoubleprime_exact':str(res2)})
    slopes=[]
    with localcontext() as ctx:
        ctx.prec=60
        for (a,ra),(b,rb) in zip(ratios,ratios[1:]):
            r=rb/ra
            alpha=-(Decimal(r.numerator)/Decimal(r.denominator)).ln()/(Decimal(b)/Decimal(a)).ln()
            slopes.append({'pair':[a,b],'definition':'-log(R_b/R_a)/log(b/a)','local_decay_slope':format(alpha,'.10f')})
    out={'issue':817,'p_exact':str(p),'systems':results,'local_slopes_only':slopes,
         'exponent_fit_performed':False,'new_enumeration_performed':False,
         'raw_pair_histograms_reconstructed':False,'L6_L7_included':False,
         'scope':'Exact revalidation of quoted atlas aggregates against independently differentiated C tables; not a rerun of the underlying enumeration.'}
    (ROOT/'result.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n');print(json.dumps(out,ensure_ascii=False,indent=2))

if __name__=='__main__':
    main()
