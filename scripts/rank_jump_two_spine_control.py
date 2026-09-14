from collections import Counter,defaultdict
import json, argparse
from pathlib import Path

DIRS=((1,0),(-1,0),(0,1),(0,-1))

def torus_rank(mask,L):
    N=L*L
    def idx(x,y): return (y%L)*L+(x%L)
    seen=set(); basis=[]
    def add(v):
        nonlocal basis
        if v==(0,0): return
        if not basis: basis=[v]; return
        a=basis[0]
        if a[0]*v[1]-a[1]*v[0]!=0: basis=[a,v]
    for root in range(N):
        if not ((mask>>root)&1) or root in seen: continue
        lift={root:(0,0)}; stack=[root]; seen.add(root)
        while stack:
            u=stack.pop(); x=u%L; y=u//L; ux,uy=lift[u]
            for dx,dy in DIRS:
                v=idx(x+dx,y+dy)
                if not ((mask>>v)&1): continue
                prop=(ux+dx,uy+dy)
                if v not in lift:
                    lift[v]=prop; seen.add(v); stack.append(v)
                else:
                    dd=(prop[0]-lift[v][0],prop[1]-lift[v][1])
                    if dd!=(0,0):
                        assert dd[0]%L==0 and dd[1]%L==0
                        add((dd[0]//L,dd[1]//L))
                        if len(basis)==2: return 2
    return len(basis)

def components(mask,L):
    N=L*L
    def idx(x,y): return (y%L)*L+(x%L)
    seen=set(); out=[]
    for root in range(N):
        if not ((mask>>root)&1) or root in seen: continue
        comp=0; stack=[root]; seen.add(root)
        while stack:
            u=stack.pop(); comp|=1<<u; x=u%L; y=u//L
            for dx,dy in DIRS:
                v=idx(x+dx,y+dy)
                if (mask>>v)&1 and v not in seen:
                    seen.add(v); stack.append(v)
        out.append(comp)
    return out

def classify(L):
    N=L*L
    neigh=sorted({(dy%L)*L+(dx%L) for dx,dy in DIRS})
    stats=Counter(); examples={}; occupancy=defaultdict(Counter)
    for rest in range(1<<(N-1)):
        m=0
        for i in range(1,N):
            if (rest>>(i-1))&1: m|=1<<i
        r0=torus_rank(m,L); r1=torus_rank(m|1,L)
        if r1-r0!=2: continue
        assert r0==0 and r1==2
        touching=[]
        for comp in components(m,L):
            att=[u for u in neigh if (comp>>u)&1]
            if att: touching.append((comp,att,torus_rank(comp|1,L)))
        single_rank2=[x for x in touching if x[2]==2]
        if single_rank2:
            mina=min(len(x[1]) for x in single_rank2)
            typ='T3' if mina==3 else 'T4plus'
        else:
            rank1=[x for x in touching if x[2]==1]
            typ='Rsplit' if len(rank1)>=2 else 'OTHER'
        stats[typ]+=1
        occupancy[typ][m.bit_count()]+=1
        examples.setdefault(typ,{
            'closed_mask':m,'open_mask':m|1,
            'touching':[{'mask':c,'attachments':att,'rank_with_v':rr} for c,att,rr in touching]
        })
    return {
        'L':L,'neighbors_of_v0':neigh,'jump2_total':sum(stats.values()),
        'type_counts':dict(stats),
        'occupancy_histograms':{t:dict(sorted(h.items())) for t,h in occupancy.items()},
        'examples':examples
    }

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--output',type=Path,required=True); a=ap.parse_args()
    out={
      'scope':'exact finite single-site rank-jump-two attachment-spine control; counts are configuration counts, not arm exponents',
      'classification':{
        'T3':'one closed outside component touches v through exactly three NN attachment germs and that component plus v already has ambient rank two',
        'T4plus':'one outside component with >=4 contacts plus v already has rank two',
        'Rsplit':'no single touching component has rank two, but at least two touching components individually create rank-one cycles through v'
      },
      'systems':[classify(3),classify(4)]
    }
    a.output.write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps(out,indent=2))
