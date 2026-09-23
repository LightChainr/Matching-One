#!/usr/bin/env python3
"""Exact parameter classification of passive width-four rank monitoring.

Consumes the unchanged PR708 certificate. Python 3.10+ and SymPy; tested with
Python 3.13.5 / SymPy 1.14.0. No sampling, new width, network access, or CI claim.
The mapping/oracle helpers are reused from the preceding passive-monitoring
analysis, not represented as a new independently authored rank implementation.

Run: python check.py input/width4-rank-closure-certificate.json --out result.json
The full default run recomputes two symbolic determinants and their coverage.
"""
from __future__ import annotations
import argparse, json, time
from collections import Counter, deque
from functools import lru_cache
from hashlib import sha1
from itertools import product
from math import isqrt
from pathlib import Path
import sympy as sp
from sympy.polys.rings import ring
from sympy.polys.domains import ZZ, QQ
from sympy.polys.matrices import DomainMatrix

R, X = ring('x', ZZ)
x, p = sp.symbols('x p')

def require(test, message):
    if not test:
        raise AssertionError(message)


def labels(keys):
    index, out = {}, []
    for key in keys:
        key = tuple(key)
        if key not in index:
            index[key] = len(index)
        out.append(index[key])
    return out


def state_for(word, transitions, initial):
    s = initial[word[0]]
    for b in word[1:]:
        s = transitions[s][b]
    return s


def physical_histories(q, initial):
    words = {}
    for a, b in product(range(16), repeat=2):
        words.setdefault(q[initial[a]][b], [a, b])
    todo = deque(words)
    while todo:
        s = todo.popleft()
        for b, t in enumerate(q[s]):
            if t not in words:
                words[t] = words[s] + [b]
                todo.append(t)
    require(len(words) == len(q), 'physical reachability')
    return [words[s] for s in range(len(q))]


def d4_quotient(q, output, initial):
    histories = physical_histories(q, initial)
    actions = []
    for sign, shift in product((1, -1), range(4)):
        masks = [sum(((b >> x) & 1) << ((sign*x+shift) % 4)
                     for x in range(4)) for b in range(16)]
        action = [state_for([masks[b] for b in word], q, initial)
                  for word in histories]
        require(len(set(action)) == len(q), 'D4 not a permutation')
        for s in range(len(q)):
            require(output[action[s]] == output[s], 'D4 output')
            for b in range(16):
                require(action[q[s][b]] == q[action[s]][masks[b]], 'D4 transition')
        actions.append(action)
    lab = labels((min(a[s] for a in actions),) for s in range(len(q)))
    representatives = [lab.index(i) for i in range(len(set(lab)))]
    weighted = [Counter((b.bit_count(), lab[q[s][b]]) for b in range(16))
                for s in representatives]
    for s in range(len(q)):
        require(Counter((b.bit_count(), lab[q[s][b]]) for b in range(16)) == weighted[lab[s]],
                'all-p strong quotient')
    return lab, representatives, weighted, [output[s] for s in representatives], histories


def graph_rank(rows):
    """Independent physical lifted-graph DFS; not certificate transitions."""
    w,h,positions,first = 4,len(rows),{},None
    for y0 in range(h):
        for x0 in range(w):
            root = y0*w+x0
            if not(rows[y0]>>x0&1) or root in positions:
                continue
            positions[root]=(0,0);stack=[root]
            while stack:
                u=stack.pop();x,y=u%w,u//w;px,py=positions[u]
                for dx,dy in ((1,0),(-1,0),(0,1),(0,-1)):
                    nx,ny=(x+dx)%w,(y+dy)%h
                    if not(rows[ny]>>nx&1):
                        continue
                    v=ny*w+nx;proposed=(px+dx,py+dy)
                    if v not in positions:
                        positions[v]=proposed;stack.append(v)
                    else:
                        a,b=proposed[0]-positions[v][0],proposed[1]-positions[v][1]
                        require(a%w==0 and b%h==0, 'physical displacement')
                        if a or b:
                            if first is None:
                                first=(a,b)
                            elif first[0]*b-first[1]*a:
                                return 2
    return int(first is not None)


def poly(cs):
    return sum((int(c)*X**k for k,c in enumerate(cs)), R.zero)


def coeff(e):
    return [int(e.get((k,), 0)) for k in range(e.degree()+1)] if e else []


def prime_check(m):
    return m >= 2 and all(m % k for k in range(2,isqrt(m)+1))


def eval_mod(cs, t, m):
    a = 0
    for c in reversed(cs):
        a = (a*t+int(c)) % m
    return a


def determinant_mod(a, m):
    a = [[int(v)%m for v in row] for row in a]
    d = 1
    for k in range(len(a)):
        i = next((i for i in range(k,len(a)) if a[i][k]), None)
        if i is None:
            return 0
        if i != k:
            a[i],a[k] = a[k],a[i]; d = -d
        pivot = a[k][k]; d = d*pivot % m; inv = pow(pivot,-1,m)
        for i in range(k+1,len(a)):
            c = a[i][k]*inv % m
            if c:
                for j in range(k+1,len(a)):
                    a[i][j] = (a[i][j]-c*a[k][j]) % m
            a[i][k] = 0
    return d % m


def rank_mod(a, m):
    a = [[int(v)%m for v in row] for row in a]; k = 0
    for j in range(len(a[0])):
        i = next((i for i in range(k,len(a)) if a[i][j]), None)
        if i is None:
            continue
        a[i],a[k] = a[k],a[i]; inv = pow(a[k][j],-1,m)
        for i in range(k+1,len(a)):
            c = a[i][j]*inv % m
            if c:
                a[i] = [(v-c*w)%m for v,w in zip(a[i],a[k])]
        k += 1
        if k == len(a):
            break
    return k


def observer_mod(weighted, output, m, t):
    """Exact finite-field lower bound; return a nonzero actual word minor."""
    n = len(output)
    rows = [[(j,mult*pow(t,k,m)%m) for (k,j),mult in row.items()] for row in weighted]
    echelon, pivots, vectors, words = [], [], [], []
    def add(v, word):
        u = v[:]
        for i,b in zip(pivots,echelon):
            c = u[i]
            if c:
                u = [(z-c*y)%m for z,y in zip(u,b)]
        i = next((i for i,z in enumerate(u) if z), None)
        if i is None:
            return
        inv = pow(u[i],-1,m)
        echelon.append([z*inv%m for z in u]); pivots.append(i)
        vectors.append(v); words.append(word)
    for b in sorted(set(output)):
        add([int(z==b) for z in output], [b])
    k = 0
    while k < len(vectors) and len(vectors) < n:
        av = [sum(c*vectors[k][j] for j,c in row)%m for row in rows]
        for b in sorted(set(output)):
            add([z if output[i]==b else 0 for i,z in enumerate(av)], [b]+words[k])
        k += 1
    minor = [[v[i] for v in vectors] for i in pivots]
    d = determinant_mod(minor,m)
    require(d != 0, 'word-minor lower bound failed')
    return {'dimension':len(vectors), 'prime':m, 'root':t, 'minor_determinant_mod':d,
            'maximum_word_length':max(map(len,words)), 'minor_rows':pivots,
            'words':words}


def determinant_poly(matrix):
    """Fraction-free Bareiss, exact divisions, full low-degree pivoting."""
    a = [row[:] for row in matrix]; n = len(a); previous = R.one; sign = 1
    if n == 1:
        return a[0][0]
    for k in range(n-1):
        choices = [(i,j) for i in range(k,n) for j in range(k,n) if a[i][j]]
        if not choices:
            return R.zero
        i,j = min(choices, key=lambda ij:(a[ij[0]][ij[1]].degree(),len(a[ij[0]][ij[1]]),ij))
        if i != k:
            a[k],a[i] = a[i],a[k]; sign = -sign
        if j != k:
            for row in a:
                row[k],row[j] = row[j],row[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                v = pivot*a[i][j]-a[i][k]*a[k][j]
                a[i][j] = v.exquo(previous) if v else v
            a[i][k] = R.zero
        previous = pivot
    return sign*a[-1][-1]


def block_order_check(matrix, blocks):
    n = len(matrix)
    require(sorted(i for b in blocks for i in b)==list(range(n)), 'not a partition')
    label = {i:k for k,b in enumerate(blocks) for i in b}
    graph = [set() for _ in blocks]
    for i,row in enumerate(matrix):
        for j,e in enumerate(row):
            if e and label[i]!=label[j]:
                graph[label[i]].add(label[j])
    indeg = [0]*len(blocks)
    for row in graph:
        for j in row:
            indeg[j] += 1
    todo = deque(i for i,v in enumerate(indeg) if not v); visited = 0
    while todo:
        i = todo.popleft(); visited += 1
        for j in graph[i]:
            indeg[j] -= 1
            if not indeg[j]: todo.append(j)
    require(visited==len(blocks), 'block graph has a cycle')


def remove_x_factors(a):
    """Only nonzero factors on x>0 are removed, never a conjectured factor."""
    a = [r[:] for r in a]; removed = 0; changed = True
    while changed:
        changed = False
        for transposed in (False,True):
            if transposed: a = list(map(list,zip(*a)))
            for i,row in enumerate(a):
                val = min(min(k[0] for k in e) for e in row if e)
                if val:
                    a[i] = [e.exquo(X**val) if e else e for e in row]
                    removed += val; changed = True
            if transposed: a = list(map(list,zip(*a)))
    return a,removed


def irreducibility_certificate(f, records):
    """Intersect possible rational factor degrees across good finite reductions."""
    possible = set(range(1,f.degree()))
    for rec in records:
        m = rec['prime']; require(prime_check(m) and int(f.LC())%m, 'bad reduction')
        fm = sp.Poly(f.as_expr(),x,modulus=m)
        require(fm.gcd(fm.diff()).degree()==0, 'non-squarefree reduction')
        _, factors = fm.factor_list()
        degrees = sorted(g.degree() for g,e in factors for _ in range(e))
        require(degrees==sorted(rec['degrees']), 'modular factor degrees disagree')
        sums = {0}
        for d in degrees: sums |= {v+d for v in sums}
        possible &= sums
    require(not possible, 'rational irreducibility not certified')


def algebraic_field(cs):
    f = sp.Poly.from_list(cs[::-1],x,domain=ZZ)
    field = QQ.algebraic_field(sp.CRootOf(f.as_expr(),0))
    return field, field.from_sympy(field.ext.as_expr())


def field_value(cs, t, field):
    result = field.zero
    for c in reversed(cs): result = result*t+field.convert(QQ(c))
    return result


def invisible_check(rec, weighted, ranks):
    F,t = algebraic_field(rec['factor']); n = len(ranks)
    v = [F.zero]*n
    for i,cs in zip(rec['support'],rec['vectors']): v[i] = field_value(cs,t,F)
    lam = field_value(rec['eigenvalue'],t,F)
    require(v[rec['support'][-1]]==F.one and sum(v,F.zero)==F.zero, 'bad contrast')
    require(all(ranks[i]!=1 for i in rec['support']), 'contrast not in binary-zero sector')
    out = [F.zero]*n
    for i,row in enumerate(weighted):
        for (k,j),m in row.items(): out[j] += v[i]*m*t**k
    require(all(a==lam*b for a,b in zip(out,v)), 'invisible eigenmode identity failed')
    zero_mass = sum((v[i] for i in rec['support'] if ranks[i]==0),F.zero)
    if rec['degree']==11:
        require(all(ranks[i]==2 for i in rec['support']), 'ternary invisibility changed')
    else:
        require(zero_mass!=F.zero, 'full-rank readout does not see claimed contrast')
    # Freeze the preparation contrast; perturb only the next row probability.
    derivative=F.zero
    for i,row in enumerate(weighted):
        for (k,j),m in row.items():
            if ranks[j]!=1:
                term=(k*t**(k-1)/(1+t)**2 if k else F.zero)-4*t**k/(1+t)**3
                derivative+=v[i]*m*term
    require(derivative!=F.zero, 'one-row detuning did not open the blind direction')
    return {'degree':rec['degree'], 'support_size':len(rec['support']),
            'one_row_detuning_first_derivative_nonzero':True,
            'eigenvector_identity_exact':True, 'ternary_invisible':rec['degree']==11,
            'histories_have_only_binary_output_zero_initially':True}


def positive_delayed(weighted, index):
    """Explicit stochastic K=UV at the cubic point and at p=1/2."""
    n = len(weighted); F,t = algebraic_field([-1,-3,-2,2]); den = (1+t)**4
    A = [[sum((m*t**k for (k,s),m in row.items() if s==j),F.zero)/den
          for j in range(n)] for row in weighted]
    z0,z6,z7,z19,z73,z74,z75 = [index(h) for h in
       ([0,0],[0,7],[0,15],[7,0],[7,0,7],[7,0,11],[7,0,13])]
    ids = [z0,z6,z7,z19,z73,z74,z75]; outside = [i for i in range(n) if i not in ids]
    require(all(j in ids for i in ids for k,j in weighted[i]), 'seven-state core not closed')
    aa=6*t*t+4*t+1; b=t**3; s=aa+3*b; c=2*t*t-2*t
    V=[]; row=[F.zero]*n; row[z0]=aa/(aa+4*b); row[z6]=4*b/(aa+4*b); V.append(row)
    row=[F.zero]*n; row[z7]=F.one; V.append(row)
    for i in (z73,z74,z75):
        row=[F.zero]*n
        for j in (z19,z73,z74,z75): row[j]=A[i][j]*den/s
        V.append(row)
    V.extend(A[i] for i in outside)
    U=[[F.zero]*92 for _ in range(n)]
    for i in (z0,z6): U[i][0]=A[i][z0]+A[i][z6]; U[i][1]=A[i][z7]
    U[z7][1]=F.one
    for k,i in enumerate((z73,z74,z75)): U[i][1]=A[i][z7]; U[i][k+2]=s/den
    U[z19][1]=A[z19][z7]
    for k,m in ((2,1),(3,2),(4,1)): U[z19][k]=m*s/(c*den)
    for k,i in enumerate(outside): U[i][5+k]=F.one
    require(all(sum(row,F.zero)==F.one for row in U+V), 'stochastic normalization')
    require(DomainMatrix(U,(n,92),F)*DomainMatrix(V,(92,n),F)==DomainMatrix(A,(n,n),F), 'K=UV')
    # Every nonzero displayed weight is positive at the unique cubic root t>1.
    cubic=sp.Poly(2*x**3-2*x*x-3*x-1,x)
    require(cubic.count_roots(1,sp.oo)==1 and cubic.count_roots(0,1)==0, 'positivity interval')
    # At half, remove the middle future-law row via an exact convex relation.
    Ah=[[QQ(sum(m for (k,s),m in row.items() if s==j),16) for j in range(n)] for row in weighted]
    require(all(Ah[z74][j]==(Ah[z73][j]+Ah[z75][j])/2 for j in range(n)), 'half convex identity')
    return {'cubic_positive_edge_states':92,'half_positive_edge_states':93,
            'cubic_factorization_shapes':[[94,92],[92,94]],
            'cubic_UV_identity_exact':True,'seven_state_core':ids,
            'positivity_proof':'all formulas positive for the unique cubic odds root t>1'}


def root_data(f):
    degree=f.degree(); cs=[int(v) for v in f.all_coeffs()][::-1]
    real=[r for r in f.intervals(eps=QQ(1,10**18)) if r[0][0]>0]
    require(len(real)==1 and real[0][1]==1, 'exception does not have one positive simple root')
    lo,hi=real[0][0]
    fp=sp.Poly(sum(cs[k]*p**k*(1-p)**(degree-k) for k in range(degree+1)),p)
    return {'odds_polynomial':cs,'p_polynomial':[int(v) for v in fp.all_coeffs()][::-1],
            'p_interval':[str(lo/(1+lo)),str(hi/(1+hi))],
            'p_decimal_diagnostic':str(sp.N((lo/(1+lo)+hi/(1+hi))/2,20))}


def run(path, seed):
    start=time.monotonic(); raw=path.read_bytes()
    blob=sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
    require(blob==seed['input_blob'], 'wrong source certificate')
    cert=json.loads(raw);q=cert['quotient_transitions'];initial=cert['quotient_initial']
    lab,reps,weighted,ranks,histories=d4_quotient(q,cert['quotient_rank_output'],initial)
    require(len(ranks)==94, 'wrong D4 state count'); bit=[int(r==1) for r in ranks]; n=94
    index=lambda h:lab[state_for(h,q,initial)]
    rows=[[(j,m*X**k) for (k,j),m in row.items()] for row in weighted]
    @lru_cache(None)
    def value(word):
        if len(word)==1: return [R.one if b==word[0] else R.zero for b in bit]
        v=value(word[1:])
        return [sum((w*v[j] for j,w in row),R.zero) if bit[i]==word[0] else R.zero
                for i,row in enumerate(rows)]
    words=list(map(tuple,seed['primary_words'])); perm=seed['column_permutation']
    require(sorted(perm)==list(range(n)), 'column permutation')
    matrix=[[value(words[perm[j]])[i] for j in range(n)] for i in range(n)]
    block_order_check(matrix,seed['observer_blocks'])
    ir={(r['block'],r['degree']):r['records'] for r in seed['irreducibility']}
    cov={(r['block'],r['degree']):r for r in seed['coverage']}
    checked=[]; exceptional={}; determinants=[]; certificate_factors=[]
    for block,ids in enumerate(seed['observer_blocks']):
        B,removed=remove_x_factors([[matrix[i][j] for j in ids] for i in ids])
        print('exact word determinant',block,len(B),flush=True)
        d=determinant_poly(B); require(bool(d), 'primary determinant identically zero')
        D=sp.Poly(d.as_expr(),x,domain=ZZ); c,fac=D.factor_list()
        reconstructed=sp.Poly(c,x,domain=ZZ)
        for f,e in fac: reconstructed*=f**e
        require(reconstructed==D, 'factorization product')
        determinants.append({'block':block,'size':len(B),'degree':D.degree(),'removed_x_power':removed})
        certificate_factors.append({'block':block,'constant':int(c),'factors':[{'coefficients':[int(v) for v in f.all_coeffs()][::-1],'multiplicity':e} for f,e in fac]})
        for f,e in fac:
            cs=[int(v) for v in f.all_coeffs()][::-1]; degree=f.degree()
            record={'block':block,'degree':degree,'multiplicity':e}
            if all(v>=0 for v in cs):
                record['reason']='nonnegative coefficients, nonzero on positive odds'; checked.append(record); continue
            if degree<=20 and f.count_roots(0,sp.oo)==0:
                record['reason']='exact Sturm count: no positive roots'; checked.append(record); continue
            irreducibility_certificate(f,ir[block,degree])
            witness=cov[block,degree]; m,t=witness['prime'],witness['root']
            require(prime_check(m) and t not in (0,m-1) and eval_mod(cs,t,m)==0, 'invalid finite reduction')
            obs=observer_mod(weighted,bit,m,t)
            require(obs['dimension']==witness['dimension'], 'modular dimension changed')
            record.update({'dimension_lower_bound':obs['dimension'],'prime':m,'root':t,
                           'minor_determinant_mod':obs['minor_determinant_mod']})
            if obs['dimension']==93:
                inv=next(r for r in seed['invisible'] if r['degree']==degree)
                require(inv['factor']==cs, 'wrong invisible factor')
                record['upper_bound']=invisible_check(inv,weighted,ranks)
                exceptional[degree]=root_data(f)
            else: require(obs['dimension']==94, 'unresolved factor')
            checked.append(record)
    require(set(exceptional)=={11,17}, 'exception set changed')
    # Full ternary readout repairs the degree-17 contrast, not degree-11.
    rescue=seed['full_rank_readout_rescue'];rr=observer_mod(weighted,ranks,rescue['prime'],rescue['root'])
    require(rr['dimension']==94, 'ternary rescue failed')
    # Independently classified transition singularities.
    A=[[sum((m*X**k for (k,s),m in row.items() if s==j),R.zero) for j in range(n)] for row in weighted]
    block_order_check(A,seed['K_components']); dk=R.one
    for ids in seed['K_components']: dk*=determinant_poly([[A[i][j] for j in ids] for i in ids])
    require(dk==X**192*(X-1)*(X+1)**27*(2*X**3-2*X**2-3*X-1)**2, 'transition determinant')
    transition=[]
    for m,t,target in ((1009,1,93),(1009,481,92)):
        am=[[eval_mod(coeff(e),t,m) for e in row] for row in A]
        am2=[[sum(am[i][k]*am[k][j] for k in range(n))%m for j in range(n)] for i in range(n)]
        require(rank_mod(am,m)==target and rank_mod(am2,m)==target, 'zero-mode rank lower bound')
        transition.append({'prime':m,'odds_mod_prime':t,'rank_A_lower':target,'rank_A_squared_lower':target})
    delayed=positive_delayed(weighted,index)
    # A source-readout mapping check on actual lifted physical graphs, not a new census.
    physical_words={tuple(histories[reps[i]]) for rec in seed['invisible'] for i in rec['support']}
    physical_words.update(tuple(h) for h in ([0,0],[0,7],[0,15],[7,0],[7,0,7],[7,0,11],[7,0,13]))
    physical_checks=0
    for h in sorted(physical_words):
        for suffix in ((),*[(b,) for b in range(16)]):
            require(graph_rank(list(h+suffix))==ranks[index(list(h+suffix))], 'physical history mapping')
            physical_checks+=1
    # Same-height preparation, avoiding a free history-length discriminator.
    states=set(initial)
    for _ in range(4):states={q[i][b] for i in states for b in range(16)}
    require(len(states)==509 and len({lab[i] for i in states})==94, 'height-five preparation')
    summary={'schema':'matching-one.p809-parameter-observability.v1','input_blob':blob,
       'scope':'fixed width-four NN site; homogeneous future Bernoulli rows; rank at virtual periodic closure',
       'binary_observable':'1{rank=1}, with current output B0 included',
       'binary_linear_dimension':'93 at the two listed roots; 94 elsewhere in 0<p<1',
       'ternary_linear_dimension':'93 only at the degree-11 root; 94 elsewhere in 0<p<1',
       'binary_exception_roots':exceptional,'primary_determinants':determinants,
       'factor_coverage':checked,'ternary_rescue_minor':{k:v for k,v in rr.items() if k!='words'},
       'det_K_in_p':'p^192*(1-p)^150*(2*p-1)*(2*p^3+p^2-1)^2',
       'transition_rank_checks':transition,'delayed_positive_models':delayed,
       'physical_graph_rank_checks':physical_checks,
       'same_height_five_preparation_states':94,
       'positive_minimality_boundary':'At observable dimension 93, minimal positive order is only bounded 93..94 unless separately constructed; do not identify linear rank with positive order.',
       'not_run':['full repository CI','new width','Monte Carlo','external independent review','publication novelty certification'],
       'sympy_version':sp.__version__}
    print('PASS all-parameter classification; elapsed seconds',time.monotonic()-start,flush=True)
    return summary,certificate_factors


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('input',type=Path)
    ap.add_argument('--certificate',type=Path,default=Path(__file__).with_name('certificate.json'))
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--factor-output',type=Path)
    args=ap.parse_args()
    result,factors=run(args.input,json.loads(args.certificate.read_text()))
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    if args.factor_output:args.factor_output.write_text(json.dumps(factors,indent=2)+'\n')
