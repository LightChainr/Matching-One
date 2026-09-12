# Continuum thin tori are not fixed-width lattice cylinders

Date: 2026-09-13. Targeted review of PR #734; no broad novelty survey.

## 1. The specific correction

PR #734, head bb41df7da8021d922a48ad353709f35912977c37, section 1a identifies
Morin-Duchesne--Saint-Aubin's thin-torus asymptotics with a fixed microscopic width
three/four lattice cylinder. Those limits are different.

The primary paper arXiv:0812.2925v2, physical PDF page 5, before equations (3)--(7),
specifies the thermodynamic limit with mesh tending to zero at critical temperature.
The next section then sends the CONTINUUM modular parameter to the cusp. The text
and displayed formulas were inspected from a rendered screenshot as well as parsed
text. Its continuum width 1 is not one (or four) occupied lattice columns.

For mesh delta on a continuum torus of periods 1 and ir, the microscopic width is
w=1/delta, and m approximately r/delta. Taking delta->0 at fixed r sends both
microscopic periods to infinity. Fixing w=4 and sending m->infinity never takes
this inner continuum limit. The latter has rank-1 concentration at EVERY fixed
p in (0,1), whereas the cited continuum formula is a critical-temperature formula.
The broad resemblance of rare wrapping in thin shapes does not identify the rates
or authorize Kac labels for the fixed-width eigenvalues.

Keep #734's correct historical identification of the ambient homology observable;
do not turn this correction into a new-observable claim. Its statement that the
rank/direction dictionary is still unmade is also stale relative to existing
spiral/cross corrections. AET 2022 computes specified-direction wrapping; the
both-projections event includes rank-one spirals and is not automatically P2.
A single directional marginal cannot by itself supply P0,P1,P2.

## 2. A precise condition for transporting a continuum exponential

Suppose, conditionally, a fixed-width lattice survival eigenvalue has expansion

    log lambda_w = -gamma/w + a/w^(1+omega) + o(w^(-1-omega)),

with a != 0 and omega>0. For m=w*r the ratio to the leading continuum exponential is

    lambda_w^m / exp(-gamma*r)
      = exp[a*r/w^omega + o(r/w^omega)].                   (1)

Relative agreement requires r=o(w^omega), assuming the remainder is controlled
uniformly along the chosen sequence. At r approximately c*w^omega there is a
nontrivial exp(a*c) factor. At still larger aspect ratios the relative mismatch can
diverge. At each fixed r the same lattice family has a perfectly correct continuum
limit. No contradiction exists.

An exact positive control is lambda_w=exp(-gamma/w+a/w^(1+omega)) for large w.
It has no fitted exponents and satisfies (1) without a remainder. Pointwise width
convergence does not imply width-uniform control after raising to the m-th power.

## 3. Root convergence likewise need not give continuum sector balance on long strips

Assume, again CONDITIONALLY and not as a square-site theorem,

    q_w-p_c = A*w^(-Delta)(1+o(1)),
    h_w'(xi_w)=B*w^(y-1)(1+o(1)),
    h_w(p)=log(lambda_2,w(p)/lambda_0,w(p)),  h_w(q_w)=0,

where xi_w lies between p_c and q_w, A and B are nonzero, and Delta>y. Also assume
the finite-m observable prefactor ratio tends to one and the observable remainder
in the log probability ratio tends to zero along the sequence under consideration.
Then the mean-value theorem gives

    log(P2(p_c)/P0(p_c))
      = -A*B*m*w^(-Delta+y-1)(1+o(1)) + o(1).             (2)

If m approximately s*w^(Delta+1-y), the conditional rank-2 probability converges to

    1/[1+exp(A*B*s)],

not generally 1/2, although q_w->p_c. In aspect-ratio language the crossover is
r approximately w^(Delta-y).

Inserting Delta=4 and y=3/4 would give r approximately w^(13/4). This arithmetic
is only an implication of the displayed hypotheses, NOT a proof of those square-site
exponents, of an all-width spectral expansion, or of uniform prefactors. The remaining
uniformity assumptions are exactly what a literature/theory bridge must provide.
A numerical cylinder-root sequence alone cannot certify them.

## 4. Primary source audit (bounded)

1. Morin-Duchesne and Saint-Aubin, Critical exponents for the homology of
   Fortuin-Kasteleyn clusters on a torus, arXiv:0812.2925v2, PRE 80 021130 (2009).
   PRIMARY_TEXT_READ: PDF pages 2--6, especially page 5 equations (3)--(7) and the
   start of section 2; page 5 screenshot checked. This is a targeted reading, not
   a claim to have checked every argument in all 19 pages.
   https://arxiv.org/pdf/0812.2925
2. Akhunzhanov, Eserkepov and Tarasevich, Exact percolation probabilities for a
   square lattice: Site percolation on a plane, cylinder, and torus,
   arXiv:2204.01517v1. PRIMARY_TEXT_READ: abstract, introduction's event definitions,
   sections 2.1--2.2. The stated L<=12 torus reach is that PAPER'S result, not a
   verified claim of the world frontier in 2026.
   https://arxiv.org/html/2204.01517v1

No originality certificate follows from this bounded reading. Equations (1)--(2)
are elementary conditional deductions supplied here, not statements attributed to
these papers. Existing uniform-aspect/root proofs remain separate from these
conditional continuum-to-lattice rate comparisons.
