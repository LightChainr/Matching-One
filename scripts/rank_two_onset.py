#!/usr/bin/env python3
"""Small controls for the all-L two-cycle-core onset proof (no new large census).

Axis periods (L,0),(0,L); diamond periods (L,L),(L,-L).
The theorem is in the companion note. Enumeration is only an existing-size
check, not the proof. Physical integer lift traversal is independent of #708.
"""
from __future__ import annotations
import argparse
from itertools import combinations
import json
from pathlib import Path
import time


def geometry(L, diamond):
    if type(L) is not int or L<2:raise ValueError('L>=2 required')
    width=2*L if diamond else L
    def canonical(x,y):
        if diamond:x-=(y//L)*L
        return (x%width,y%L)
    sites=[(x,y) for y in range(L) for x in range(width)]
    index={xy:i for i,xy in enumerate(sites)}
    adj=[]
    for x,y in sites:
        adj.append([(index[canonical(x+dx,y+dy)],dx,dy)
                    for dx,dy in ((1,0),(-1,0),(0,1),(0,-1))])
    return sites,index,canonical,adj


def rank(mask,adj):
    positions={};first=None
    for root in range(len(adj)):
        if not (mask>>root)&1 or root in positions:continue
        positions[root]=(0,0);todo=[root]
        while todo:
            i=todo.pop();x,y=positions[i]
            for j,dx,dy in adj[i]:
                if not (mask>>j)&1:continue
                proposed=(x+dx,y+dy)
                if j not in positions:positions[j]=proposed;todo.append(j)
                else:
                    v=(proposed[0]-positions[j][0],proposed[1]-positions[j][1])
                    if v!=(0,0):
                        if first is None:first=v
                        elif first[0]*v[1]-first[1]*v[0]:return 2
    return int(first is not None)


def predicted_minimizers(L,diamond):
    sites,index,canonical,_=geometry(L,diamond);out=set()
    def encode(points):return sum(1<<index[p] for p in set(points))
    if not diamond:
        for x in range(L):
            for y in range(L):
                out.add(encode([(j,y) for j in range(L)]+[(x,j) for j in range(L)]))
    else:
        # Full straight physical row/column plus either length-L transverse arc.
        # Set deduplication is checked against the theorem's 4 L^2 count.
        for x,y in sites:
            for direction in (0,1):
                full=[canonical(x+j,y) if direction==0 else canonical(x,y+j) for j in range(2*L)]
                for sign in (-1,1):
                    plug=[canonical(x,y+sign*j) if direction==0 else canonical(x+sign*j,y) for j in range(1,L)]
                    out.add(encode(full+plug))
    return out


def check_small(L,diamond):
    sites,_,_,adj=geometry(L,diamond);onset=(3*L-1 if diamond else 2*L-1)
    actual=set();checked=0;below=0
    for k in range(onset+1):
        for chosen in combinations(range(len(sites)),k):
            mask=sum(1<<j for j in chosen);checked+=1
            if rank(mask,adj)==2:
                if k<onset:below+=1
                else:actual.add(mask)
    expected=predicted_minimizers(L,diamond)
    if below or actual!=expected:raise AssertionError('geometric classification failed')
    return {'L':L,'geometry':'diamond' if diamond else 'axis','sites':len(sites),
            'onset':onset,'configurations_checked':checked,'rank2_below_onset':below,
            'minimizers':len(actual),'all_minimizers_match_classification':True,
            'minimum_configuration_masks':sorted(actual)}


def report():
    started=time.perf_counter()
    cases=[check_small(2,False),check_small(3,False),check_small(4,False),
           check_small(2,True),check_small(3,True)]
    constructions=[]
    for L in range(2,9):
        for diamond in (False,True):
            states=predicted_minimizers(L,diamond);adj=geometry(L,diamond)[3]
            mass=3*L-1 if diamond else 2*L-1;expected=4*L*L if diamond else L*L
            if len(states)!=expected or any(s.bit_count()!=mass or rank(s,adj)!=2 for s in states):
                raise AssertionError('all-L construction control failed')
            constructions.append({'L':L,'geometry':'diamond' if diamond else 'axis',
                                  'constructed_configurations':len(states),'mass':mass})
    return {'schema':'matching-one.two-cycle-core-onsets.v1',
      'theorems':{'axis':{'minimum_mass':'2*L-1','minimizer_count':'L^2','structure':'one full physical row plus one full physical column'},
                  'diamond':{'minimum_mass':'3*L-1','minimizer_count':'4*L^2','structure':'one full physical row/column of 2L vertices plus L-1 internal vertices of a straight transverse length-L arc'}},
      'scope':'NN square site, honest axis/diamond quotients, L>=2; rank2 only',
      'small_existing_size_checks':cases,'total_configurations_checked':sum(c['configurations_checked'] for c in cases),
      'construction_checks':constructions,'diamond_L5_onset_prediction':{'mass':14,'count':100,'enumerated':False},
      'proof_not_enumeration':'two independent cycles reduce to a wedge or theta; systolic bounds and equality cases',
      'seconds':time.perf_counter()-started,
      'limits':['no new diamond L4/L5 enumeration','does not classify every higher-mass wrapping cell','no threshold exponent or continuum field claim']}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path);args=ap.parse_args()
    text=json.dumps(report(),indent=2,allow_nan=False)+'\n'
    if args.out:
        with args.out.open('x',encoding='utf-8') as f:f.write(text)
    else:print(text,end='')
