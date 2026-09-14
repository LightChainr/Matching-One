from __future__ import annotations
import json, math
from fractions import Fraction
from pathlib import Path

SYSTEMS = {
  'square_L3': {
    'L':3,
    'rank':[ [1,9,36,78,90,45,0,0,0,0], [0,0,0,6,36,72,48,0,0,0], [0,0,0,0,0,9,36,36,9,1] ],
    'balance_p': '0.58651145511267563565455897660690173482430062489383929'
  },
  'square_L4': {
    'L':4,
    'rank':[ [1,16,120,560,1812,4272,7448,9440,8082,3984,792,32,0,0,0,0,0], [0,0,0,0,8,96,560,1984,4580,6368,4704,1472,160,0,0,0,0], [0,0,0,0,0,0,0,16,208,1088,2512,2864,1660,560,120,16,1] ],
    'balance_p': '0.59067211233102829689590201143951286962111713272216'
  },
  'triangular_L3': {
    'L':3,
    'rank':[ [1,9,36,75,45,0,0,0,0,0], [0,0,0,9,81,81,9,0,0,0], [0,0,0,0,0,45,75,36,9,1] ],
    'balance_p': '0.5'
  },
  'triangular_L4': {
    'L':4,
    'rank':[ [1,16,120,560,1808,4128,6304,5616,2160,304,0,0,0,0,0,0,0], [0,0,0,0,12,240,1704,5520,8550,5520,1704,240,12,0,0,0,0], [0,0,0,0,0,0,0,304,2160,5616,6304,4128,1808,560,120,16,1] ],
    'balance_p': '0.5'
  }
}

def sector_prob(coeff,p):
    N=len(coeff)-1
    return sum(c*p**k*(1-p)**(N-k) for k,c in enumerate(coeff))

def expected_permutation_gap(rank1):
    N=len(rank1)-1
    return sum(Fraction(c, math.comb(N,k)) for k,c in enumerate(rank1))

def evaluate(name, obj):
    p=float(obj['balance_p']); L=obj['L']; N=L*L
    probs=[sector_prob(c,p) for c in obj['rank']]
    s=sum(probs); probs=[x/s for x in probs]
    P0,P1,P2=probs
    b=.5*math.log(P0/P2)
    d=math.log(P1/math.sqrt(P0*P2))
    chi_re=(3*P1-1)/2
    chi_im=math.sqrt(3)*(P2-P0)/2
    tb=[1,0,-1]; td=[0,1,0]
    Eb=sum(probs[i]*tb[i] for i in range(3)); Ed=P1
    Fbb=sum(probs[i]*tb[i]**2 for i in range(3))-Eb*Eb
    Fdd=P1*(1-P1)
    Fbd=sum(probs[i]*tb[i]*td[i] for i in range(3))-Eb*Ed
    gap=expected_permutation_gap(obj['rank'][1])
    return {
      'name':name,'L':L,'N':N,'balance_p':obj['balance_p'],
      'P0_P1_P2':probs,'b':b,'d':d,'chi':[chi_re,chi_im],
      'fisher_bd':[[Fbb,Fbd],[Fbd,Fdd]],
      'expected_permutation_gap_exact':f'{gap.numerator}/{gap.denominator}',
      'expected_permutation_gap':float(gap),
      'gap_over_L_5over4':float(gap)/(L**1.25),
      'mean_continuous_p_gap':float(gap)/(N+1)
    }

def main():
    out={'scope':'exact coefficient postprocessing; balance roots from existing exact controls; no scaling fit',
         'systems':[evaluate(k,v) for k,v in SYSTEMS.items()]}
    for r in out['systems']:
        P0,P1,P2=r['P0_P1_P2']; x,y=r['chi']; d=r['d']
        inv=[(1-x)/3-y/math.sqrt(3),(1+2*x)/3,(1-x)/3+y/math.sqrt(3)]
        assert max(abs(inv[i]-[P0,P1,P2][i]) for i in range(3))<2e-12
        lhs=(1+2*x)**2
        rhs=math.exp(2*d)*((1-x)**2-3*y*y)
        assert abs(lhs-rhs)<2e-11
    Path('rank-character-persistence-controls.json').write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps(out,indent=2,ensure_ascii=False))
if __name__=='__main__': main()
