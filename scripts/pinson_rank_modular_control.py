from __future__ import annotations
import math, json
import mpmath as mp
mp.mp.dps=60

# Q=1 specialization of Arguin hep-th/0111193, Eqs. (7),(9),(10),(13).
# Common Dedekind/prefactor cancels in rank probability ratios.
G0=mp.mpf(2)/3
E0=mp.mpf(2)/3

def gcd2(m,n):
    if m==0 and n==0:return 0
    return math.gcd(abs(m),abs(n))

def rank_probs(tau_r,tau_i,cut=30):
    A=mp.mpf('0'); R1=mp.mpf('0')
    for m in range(-cut,cut+1):
        for n in range(-cut,cut+1):
            w=mp.e**(-mp.pi*G0*(m*m*tau_i*tau_i+(n-m*tau_r)**2)/tau_i)
            g=gcd2(m,n)
            A += w*((-1)**g)
            if m or n:
                R1 += w*(mp.cos(mp.pi*E0*g)-((-1)**g))
    Z0=A/2; Z=A+R1
    return Z0/Z,R1/Z,Z0/Z

def diagnostics(tau_r,tau_i,cut=30):
    p0,p1,p2=rank_probs(tau_r,tau_i,cut)
    c=p1/(2*mp.sqrt(p0*p2)); discr=p1*p1-4*p0*p2
    b=mp.mpf('.5')*mp.log(p0/p2)
    if c < 1:
        th=mp.acos(-c); zeros=[mp.mpc(b,th),mp.mpc(b,-th)]; typ='complex-centered'
    elif c == 1:
        zeros=[mp.mpc(b,mp.pi),mp.mpc(b,mp.pi)]; typ='double'
    else:
        a=mp.acosh(c); zeros=[mp.mpc(b+a,mp.pi),mp.mpc(b-a,mp.pi)]; typ='negative-u-real-split'
    entropy=-sum(x*mp.log(x) for x in (p0,p1,p2))
    chi=p0*mp.e**(-2j*mp.pi/3)+p1+p2*mp.e**(2j*mp.pi/3)
    return dict(P0=p0,P1=p1,P2=p2,c=c,discriminant=discr,zero_type=typ,zeros=zeros,
                entropy=entropy,rank_character=chi)

def a_rect(r):return diagnostics(0,mp.mpf(r),35)['P0']
rstar=mp.findroot(lambda r:a_rect(r)-mp.mpf('.25'),(mp.mpf('1.4'),mp.mpf('2.2')))
rows={}
for r in (mp.mpf('.25'),mp.mpf('.5'),mp.mpf('1'),mp.mpf('2'),mp.mpf('4'),rstar):
    d=diagnostics(0,r,35)
    rows[mp.nstr(r,30)]={
      'P0':mp.nstr(d['P0'],45),'P1':mp.nstr(d['P1'],45),'P2':mp.nstr(d['P2'],45),
      'c':mp.nstr(d['c'],45),'discriminant':mp.nstr(d['discriminant'],45),
      'zero_type':d['zero_type'],
      'zeros':[{'re':mp.nstr(z.real,35),'im':mp.nstr(z.imag,35)} for z in d['zeros']],
      'entropy':mp.nstr(d['entropy'],45),
      'rank_character':{'re':mp.nstr(d['rank_character'].real,45),'im':mp.nstr(d['rank_character'].imag,45)} }
z=mp.mpc('.31','1.27');mods={}
for label,t in [('tau',z),('Ttau',z+1),('Stau',-1/z)]:
    d=diagnostics(t.real,t.imag,35)
    mods[label]=[mp.nstr(d[k],45) for k in ('P0','P1','P2')]
out={
 'scope':'Q=1 continuum torus rank probabilities from Arguin/Pinson homology sums; numerical Gaussian sums, not a rigorous interval implementation',
 'formula':'Arguin hep-th/0111193 Eqs (7),(9),(10),(13), Q=1 => g/4=2/3,e0=2/3; common boson prefactor cancels',
 'rectangular':rows,
 'critical_zero_collision_aspect_r_gt_1':mp.nstr(rstar,55),
 'reciprocal_collision_aspect':mp.nstr(1/rstar,55),
 'collision_rank_law':['1/4','1/2','1/4'],
 'modular_checks':mods,
 'max_modular_probability_difference':mp.nstr(max(abs(mp.mpf(mods['tau'][i])-mp.mpf(mods[x][i])) for x in ('Ttau','Stau') for i in range(3)),20)
}
print(json.dumps(out,indent=2))
