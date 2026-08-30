#!/usr/bin/env python3
"""Exact cyclic quotient for the C4 self-matching checkerboard triangulation.

For primitive odd (a,b), N=a^2+b^2 is even and
Z^2/<(a,b),(-b,a)> is cyclic with label j=a*x+b*y mod N.

Square NN steps are +/-a,+/-b.  Add the checkerboard diagonal connecting the
even-parity corners of every square.  In cyclic labels this is equivalent to
adding, from every even j, the two undirected diagonal edges
j--(j+a+b) and j--(j+a-b).

The square cell with lower-left label j is split into two triangles:

  j even: (j,j+a,j+a+b), (j,j+b,j+a+b)
  j odd : (j,j+a,j+b),   (j+a,j+b,j+a+b)

The validator checks the torus triangulation incidences and the exact C4
automorphism j -> t*j, t=b*a^{-1} mod N, for which t^2=-1 mod N.
"""
from __future__ import annotations

import argparse
import math
from collections import Counter
from dataclasses import dataclass
from typing import Iterable, Tuple

Edge=Tuple[int,int]
Face=Tuple[int,int,int]

def edge(u:int,v:int)->Edge:
    if u==v: raise ValueError("loop edge")
    return (u,v) if u<v else (v,u)

def face_key(x:Iterable[int])->Face:
    a=tuple(sorted(x))
    if len(set(a))!=3: raise ValueError("degenerate triangular face")
    return a  # type: ignore[return-value]

@dataclass(frozen=True)
class Geometry:
    a:int; b:int; n:int; edges:frozenset[Edge]; faces:tuple[Face,...]; rotation_multiplier:int

def build(a:int,b:int)->Geometry:
    if a<=0 or b<=0 or math.gcd(a,b)!=1 or a%2!=1 or b%2!=1:
        raise ValueError("require positive, odd, primitive (a,b)")
    n=a*a+b*b
    if n%2: raise AssertionError("odd a,b must give even N")
    es=set(); fs=[]
    for j in range(n):
        # two forward NN edges per lattice site
        es.add(edge(j,(j+a)%n)); es.add(edge(j,(j+b)%n))
        # every even site is the even-parity corner shared by four diagonals;
        # two forward choices generate every undirected diagonal once.
        if j%2==0:
            es.add(edge(j,(j+a+b)%n)); es.add(edge(j,(j+a-b)%n))

        ja=(j+a)%n; jb=(j+b)%n; jab=(j+a+b)%n
        if j%2==0:
            fs.append(face_key((j,ja,jab)))
            fs.append(face_key((j,jb,jab)))
        else:
            fs.append(face_key((j,ja,jb)))
            fs.append(face_key((ja,jb,jab)))

    inv_a=pow(a,-1,n)
    t=(b*inv_a)%n
    if (t*t+1)%n!=0:
        raise AssertionError("rotation multiplier does not square to -1")
    return Geometry(a,b,n,frozenset(es),tuple(fs),t)

def validate(g:Geometry)->dict[str,int]:
    n=g.n
    if len(g.edges)!=3*n: raise AssertionError(f"E={len(g.edges)} != 3N")
    if len(g.faces)!=2*n: raise AssertionError(f"F={len(g.faces)} != 2N")
    if len(set(g.faces))!=len(g.faces): raise AssertionError("duplicate faces")

    incidence=Counter()
    for f in g.faces:
        efs=(edge(f[0],f[1]),edge(f[0],f[2]),edge(f[1],f[2]))
        if any(e not in g.edges for e in efs): raise AssertionError("face uses missing edge")
        incidence.update(efs)
    if set(incidence)!=set(g.edges): raise AssertionError("some graph edge belongs to no face")
    if set(incidence.values())!={2}: raise AssertionError("torus edge is not incident to exactly two faces")

    degree=[0]*n
    for u,v in g.edges: degree[u]+=1; degree[v]+=1
    ev={degree[j] for j in range(n) if j%2==0}; od={degree[j] for j in range(n) if j%2==1}
    if ev!={8} or od!={4}: raise AssertionError(f"unexpected degree sets even={ev}, odd={od}")

    # A triangulation is self-matching: every pair of vertices on every face is already an edge.
    matching_added=0
    for f in g.faces:
        for i in range(3):
            for k in range(i+1,3):
                if edge(f[i],f[k]) not in g.edges: matching_added+=1
    if matching_added: raise AssertionError("matching construction would add an edge")

    # Exact C4 automorphism on the cyclic quotient.
    t=g.rotation_multiplier
    if t%2!=1: raise AssertionError("rotation must preserve checkerboard parity")
    rotated={edge((t*u)%n,(t*v)%n) for u,v in g.edges}
    if rotated!=set(g.edges): raise AssertionError("90-degree map is not a graph automorphism")
    rotated_faces={face_key(((t*x)%n for x in f)) for f in g.faces}
    if rotated_faces!=set(g.faces): raise AssertionError("90-degree map is not a face automorphism")

    # t^2=-1 and therefore t^4=1 on vertices.
    if (t*t+1)%n or (pow(t,4,n)-1)%n: raise AssertionError("C4 algebra failed")

    return {"V":n,"E":len(g.edges),"F":len(g.faces),"even_degree":8,"odd_degree":4,"rotation_multiplier":t}

def main()->int:
    p=argparse.ArgumentParser(description=__doc__); p.add_argument("a",type=int); p.add_argument("b",type=int); args=p.parse_args()
    g=build(args.a,args.b); info=validate(g)
    print(f"(a,b)=({g.a},{g.b}) N={g.n}")
    print(" ".join(f"{k}={v}" for k,v in info.items()))
    print("self-matching triangulation and C4 automorphism: PASS")
    return 0
if __name__=="__main__": raise SystemExit(main())
