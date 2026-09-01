from pathlib import Path
import importlib.util,json,itertools,hashlib
from fractions import Fraction as F
oracle_path=Path(__file__).resolve().parent.parent/'p337-face-kernel-20260831/verify.py'
spec=importlib.util.spec_from_file_location('face_oracle',oracle_path)
o=importlib.util.module_from_spec(spec);spec.loader.exec_module(o)

def support(vertices,edges,As):
 adj={v:[] for v in vertices}
 for i,(v,w,g) in enumerate(edges):adj[v].append((w,i));adj[w].append((v,i))
 tin={}; low={};bridges=set()
 def visit(v,parent_edge=-1):
  tin[v]=low[v]=len(tin)
  for w,i in adj[v]:
   if i==parent_edge:continue
   if w not in tin:
    visit(w,i);low[v]=min(low[v],low[w])
    if low[w]>tin[v]:bridges.add(i)
   else:low[v]=min(low[v],tin[w])
 for v in vertices:
  if v not in tin:visit(v)
 cyclic={v for v in As if any(i not in bridges for w,i in adj[v])}
 core=set(vertices);changed=True
 while changed:
  remove={v for v in core if sum(w in core for w,i in adj[v])<2}
  changed=bool(remove);core-=remove
 return cyclic,core,bridges

# Two abstract four-cycles joined by a two-edge bridge path through an A switch.
As={'al0','al1','ar0','ar1','bridge'};Bs={'bl0','bl1','br0','br1'}
edges=[('al0','bl0',(5,5)),('bl0','al1',(0,0)),('al1','bl1',(0,0)),('bl1','al0',(0,0)),('ar0','br0',(-5,5)),('br0','ar1',(0,0)),('ar1','br1',(0,0)),('br1','ar0',(0,0)),('bl1','bridge',(0,0)),('bridge','br0',(0,0))]
cyc,core,bridges=support(As|Bs,edges,As)
assert cyc==As-{'bridge'} and 'bridge' in core
comparisons=0
for bits in itertools.product((0,1),repeat=4):
 retained={a for a,z in zip(sorted(cyc),bits) if z}
 pair=[]
 for z in (0,1):
  occupied=Bs|retained|({'bridge'} if z else set())
  pair.append(o.graph_stats(occupied,[e for e in edges if e[0] in occupied and e[1] in occupied]))
 assert (pair[0]['beta'],pair[0]['rank'])==(pair[1]['beta'],pair[1]['rank'])
 comparisons+=1

# A fixed honest N50 B configuration: two local 4-cycles joined by one A bridge.
Bchain={o.key(0,-1),o.key(1,0),o.key(1,2),o.key(2,3)}
potential=o.A|Bchain;pes=o.edge_set(potential,o.NN_HALF)
cyclic,core,bridges=support(potential,pes,o.A)
a_bridge=o.key(1,1)
expected_cyclic={o.key(0,0),o.key(1,-1),o.key(1,3),o.key(2,2)}
assert cyclic==expected_cyclic and a_bridge in core and a_bridge not in cyclic
left=o.key(0,0)
configs=[('bridge_all_A',o.A),('bridge_removed',o.A-{a_bridge}),('one_cycle_A_removed',o.A-{left}),('one_cycle_A_and_bridge_removed',o.A-{left,a_bridge}),('all_offcycle_A_integrated_support',cyclic)]
checked=[o.check(name,aa,Bchain) for name,aa in configs]
assert [(x['parent']['beta'],x['parent']['rank'],x['S_direct']) for x in checked]==[(2,0,18),(2,0,21),(1,0,19),(1,0,22),(2,0,81)]

# Potential full-A square graph has no cyclic A articulation. Two holes can
# create a shared-A figure-eight in a selected subgraph, so dynamic splitting is unsafe.
Bcross={o.key(1,0),o.key(0,1),o.key(-1,0),o.key(0,-1)}
centre=o.key(0,0);holes={o.key(-1,1),o.key(1,-1)}
actual=o.A-holes
bow=[o.check('shared_A_on',actual,Bcross),o.check('shared_A_off',actual-{centre},Bcross)]
assert [(x['parent']['beta'],x['S_direct']) for x in bow]==[(2,24),(0,23)]

alpha=F(2,5);r=1-alpha;w=F(1,8);factor=r+alpha*w;rho=alpha*w/factor
assert factor==F(13,20) and rho==F(1,13)
block=r*r+2*r*alpha*w+alpha*alpha*w*w*4
assert block==F(43,100)
K=4+21*rho+2*(2*r*alpha*w+2*alpha*alpha*w*w*4)/block
S=89-3*21*rho+2*(-3*2*r*alpha*w-4*alpha*alpha*w*w*4)/block
beta=2*alpha*alpha*w*w*4/block
assert S==2*beta-3*K+101
out={'status':'PASS','oracle_path':str(oracle_path),'oracle_sha256':hashlib.sha256(oracle_path.read_bytes()).hexdigest(),'abstract_bridge_pairs':comparisons,'abstract_bridge_in_2core':True,'fixed_N50_B_configurations':2,'real_graph_checks':checked+bow,'bridge_cycle_A_parent_coordinates':sorted(map(o.coordinate,cyclic)),'bridge_offcycle_A_count':len(o.A-cyclic),'closed_message_example':{'alpha':str(alpha),'exp_t':'2','F':str(factor),'rho':str(rho),'two_A_block':str(block),'mean_K':str(K),'mean_S':str(S),'mean_beta':str(beta),'partition_expression':'2^89*(13/20)^21*(43/100)^2'},'new_B_scan':False,'new_random_samples':0,'cloud_jobs':0}
Path(__file__).with_name('cycle-check.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ['status','abstract_bridge_pairs','fixed_N50_B_configurations','bridge_offcycle_A_count','closed_message_example']}))
