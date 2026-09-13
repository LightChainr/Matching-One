# P398: a baseline-invisible odd sector is fully visible to two separated insertions

Date: 2026-09-12. Completed bounded exact analysis at widths 4 and 5.
Related existing channels: #594, #598, #600, #601, #610, #644. No new production,
new width campaign or generic literature ticket is requested. This is P398's
calibration process, **not square-site percolation**. No mathematical novelty
claim is made for parity, response expansions or bilinear realization theory.

## 1. The question changed; the original readout dictionary did not

The archived `notes/p398-reflection-parity-20260906.md` already establishes the
reflection-even baseline quotient, vanishing first-order response to an odd
perturbation, and an observed quadratic perturbation-size ladder. We do not
claim those again. The remaining question answered here is different:

> Which dynamics hidden at baseline can be recovered from the *delay dependence*
> of two odd insertions, with the original sources and readouts unchanged?

The distinction is experimentally material: measuring one small quadratic
coefficient is not measuring the dynamics between two interventions.

We read the implementation at main commit
`eb89e9422791d9e3c3a78f0e65d56912b815a7bd`:

- `scripts/p398_intervention_transport.py`, blob
  `ac13194751b8ba727a45499900e3f12ce1b90461`: generator, original sources/readouts;
- `scripts/planar_state_operations.py`, blob
  `c5f57dcb8606dbb62306cde839edefae680f863d`: cyclic join and point detach;
- `scripts/p398_reflection_parity.py` and the archived parity note: normalization
  and reflection convention.

The new standalone standard-library script independently implements these rules.
It is not a byte-for-byte import of the entire original module, nor a full
repository checkout/test execution.

## 2. Fixed finite model and physical perturbation

States are noncrossing partitions of w cyclically labelled points, encoded as
canonical restricted growth strings. For each point j let J_j be the row-oriented
join-move generator (merge j and j+1 modulo w), and D_j the detach generator.
Each is a deterministic move minus the identity; a no-op cancels exactly.

    G = sum_j (J_j + D_j),
    R : j -> w-1-j,
    H = J_0 - J_(w-2).

Thus RGR=G and RHR=-H. The intervention is physical: G+epsilon H has join rates
1+epsilon and 1-epsilon at the two named sites, and rate one everywhere else.
It is a Markov generator for |epsilon|<=1, and is irreducible for |epsilon|<1
at the two widths checked. H by itself is **not** a Markov generator.

Our H is twice the archived `H_odd=(J_0-J_(w-2))/2`. Every two-insertion kernel
below is therefore four times its value in the archived normalization. Ranks
and pole locations are unchanged by this nonzero rescaling.

F has the ORIGINAL three readout columns: block count, singleton count,
`wrap=1[state[0]==state[-1]]`. S has the ORIGINAL four probability-source rows:
delta all singletons, delta single block, delta wrapped pair, and the uniform
law. All are R-even. Exact computations store n*S, with the common denominator
n recorded, to avoid introducing a floating-point uniform source.

## 3. General algebra: forbidden once, allowed twice

Let P_+=(I+R)/2 and P_-=(I-R)/2. Since G commutes with R, it restricts to G_+
and G_-. An odd H maps each parity to the other. Hence for all nonnegative
waiting times,

    S exp(t0 G) H exp(t1 G) F = 0.

Inserting R on both sides proves this identity pointwise. More generally any
word with an odd number of H factors has zero even-to-even matrix element.
"Not forbidden" does NOT imply nonzero: parity neither guarantees that an
allowed response fires nor fixes its sign.

Now define the matrix-valued delay kernel

    K(tau) = S H exp(tau G) H F
           = C_- exp(tau G_-) B_-,

where B_- is H F in odd coordinates and C_- is S H restricted to them. Its
Taylor coefficients are S H G^k H F. An exact linear realization of this delay
kernel has order equal to its reachable-and-observable odd subspace, NOT
necessarily the entire odd dimension. The equality is tested below rather
than assumed from the selection rule.

### Relation to projected memory

In parity coordinates write H_+-: E_- -> E_+ and H_-+: E_+ -> E_-. For a fixed
amplitude epsilon, wherever the resolvents exist, the exact Schur complement is

    P_+ (zI-G-epsilon H)^(-1) P_+
      = [zI-G_+ - epsilon^2 H_+-(zI-G_-)^(-1)H_-+]^(-1)

on the even subspace. Eliminating the odd component from
x'=(G+u(t)H)x with x_-(0)=0 equivalently gives

    x_+'(t) = G_+ x_+(t)
      + u(t) H_+- int_0^t exp((t-s)G_-) u(s) H_-+ x_+(s) ds.

Thus the baseline quotient is exact when u=0, but the controlled even dynamics
has an odd-propagation memory kernel. No time-local even-only replacement is
asserted. This connects the archived parity and projected-memory questions
without conflating their particular numerical model reductions.

## 4. Two actual rate pulses, not a signed fictitious propagator

For duration delta>0 use two separately tunable physical windows:

    Y_delta(e1,e2;tau)
      = S exp(delta(G+e1 H)) exp(tau G) exp(delta(G+e2 H)) F.

All three factors are Markov semigroups for |e1|,|e2|<=1. Define

    J_delta = int_0^delta exp(sG) H exp((delta-s)G) ds.

Then the mixed derivative at zero is exactly

    d_e1 d_e2 Y_delta(0,0;tau)
      = S J_delta exp(tau G) J_delta F,

and J_delta/delta -> H as delta->0. Only after division by delta^2 and this
short-duration limit is the result K(tau). Finite-duration corrections must
not be ignored or relabelled as sampling error.

The four-sign contrast

    [Y(e,e)-Y(e,-e)-Y(-e,e)+Y(-e,-e)]/(4e^2)

converges to the finite-window mixed derivative with O(e^2) error. Reflection
implies Y(e,e)=Y(-e,-e) and Y(e,-e)=Y(-e,e). These equalities do not require
that G be reversible, or that the source be its stationary law.

**Finite-duration persistence.** J_delta is odd and analytic in delta. If a
d-by-d Hankel minor of K is nonzero, the same minor for the finite-window
kernel equals delta^(2d) times that nonzero minor plus O(delta^(2d+1)). Hence
the full d-mode order survives for every sufficiently small positive duration;
it is not confined to a physically unattainable zero-duration pulse. This is
an existence statement, not a noise budget or a numerical lower bound on the
allowed pulse duration.

## 5. Executed exact results

| width | microscopic n | unforced raw I/O order | two-insertion delay order | full G/H controlled word order | source-minus-uniform controlled order |
|---|---:|---:|---:|---:|---:|
| 4 | 14 | 10 | 4 | 14 | 13 |
| 5 | 42 | 26 | 16 | 42 | 41 |

The unforced raw orders exhaust the even dimensions. The two-insertion orders
exhaust the odd dimensions. The fourth column concerns only K(tau); it is NOT
added to the baseline rank by assumption. The fifth column has its own complete
noncommutative Hankel certificate for all words in G and H.

The order notion in the last two columns is the minimum dimension of a
**homogeneous linear/bilinear input-output realization**, allowing the declared
initial preparations. The original model supplies the upper bound n. For source
contrasts the constant function is invariant under G, killed by H, and invisible
to the contrasts, giving upper bound n-1. Explicit nonzero minors attain both
bounds. A model allowing an explicit affine/output offset need not count that
constant in the same way; hence both raw and contrasted results are reported.
No assertion about minimum nonlinear, positive, or microscopic state dimension
follows. All control derivatives are around zero inside the physical rate
interval; a formal word in H is not an independently executable negative-rate
transition.

### A single ORIGINAL scalar observer sees all four odd modes at w=4

Take the original delta_wrapped_pair source and wrap readout. Its delay kernel
has the exact Laplace transform

    Khat(z) = 2(z^2+11z+27) /
              [(z^2+10z+23)(z^2+11z+26)].

The first coefficients K^(k)(0) are

    0, 2, -20, 156, -1122, 7822, -53932, 371172, -2561154.

The 4x4 Hankel determinant from coefficient zero is exactly -16. The denominator
is the odd generator's characteristic polynomial, checked by integer
Faddeev-LeVerrier and exact annihilation. The numerator follows from these
initial moments. The nonzero Hankel determinant proves that no lower-order
constant linear realization reproduces the whole delay curve.

In particular K(0)=0 while K'(0)=2: a coincident-pulse measurement can be zero
although the separated-pulse experiment contains four dynamic modes. This is
why delay is an information-bearing coordinate, not a plot embellishment.
At w=5 the original all-singletons source and wrap readout alone attain order
16, certified by a nonzero 16x16 integer Hankel determinant modulo the stated
prime. This does not claim that every one of the source/readout pairs does so.
The single-block source has no initial join response and is a built-in zero.

### No universal sign of the quadratic effect

Already at w=4, all-singletons -> blocks has K(0)=2. All-singletons -> wrap has
K(0)=0 and K'(0)=-2, whereas wrapped-pair -> wrap has K(0)=0 and K'(0)=+2.
Analyticity gives opposite signs at sufficiently small positive delay for the
last two. Thus reflection parity fixes allowed orders, not the sign of a
quadratic response. These are exact examples in the existing dictionary.

## 6. Why the rank certificate is exact

The output contains integer G,H, the reflected state permutation, the original
scaled sources/readouts, exact odd/even coordinates, and explicit Hankel minors.
For efficiency, candidate independent rows and columns are selected modulo
1000000007. Every printed minor is also stored as an INTEGER matrix and its
determinant modulo this prime is nonzero. Consequently it is nonzero over Q.

Modular stabilization alone would NOT prove a rational upper bound. Here the
independent mathematical bounds are parity dimension, the microscopic n, and
the invisible constant. The script refuses a conclusion when a lower bound
does not attain the stated bound. No floating rank tolerance is used.

All ten stored rank-witness minors were also independently recomputed by
fraction-free Bareiss elimination over the integers, with exact division at
every step. Every determinant is nonzero and its modular residue agrees. The
separate verifier reconstructs the selected word coefficients before checking
the minors. The w=4 scalar determinant is additionally checked by Fraction
elimination (-16). For each odd generator, integer characteristic
polynomial coefficients divide exactly at every Faddeev-LeVerrier step, the
terminal matrix is zero, and the scalar moment recurrence is checked. Seven
small mathematical tests execute these facts, including the quarter-strength
normalization under replacement H -> H/2 and the actual positive/negative
examples. The tests are not a claim of full repository CI.

## 7. Executed physical-window numerical control

`p398_double_pulse_numerical_control.py` uses mpmath at 45 decimal digits on the
w=4 scalar observer above. At tau=1/4 and delta=1/32 it compares actual matrix
exponentials against an independent 35-term Frechet-derivative power series.
The latter is a numerical convergence control, not a rigorous remainder bound.

Finite-window mixed derivative:

    0.00011718197680195303035423582807185137.

Centered four-sign errors at e=1/8,1/16,1/32 are respectively

    6.26516675e-10, 1.56628921e-10, 3.91572148e-11,

approximately quartering when the amplitude halves, as the expansion predicts.
Reflection-related matrix-element pairs agree at the recorded precision.

At this finite duration the mixed derivative divided by delta^2 is
0.1199943442, while the ideal short-pulse K(1/4) is 0.1476852668. This visible
difference is intentional: the executable protocol and its limiting formula
are distinct objects. Exact observability likewise is not a claim that noisy
finite-amplitude data can estimate all modes cheaply. No sample top-up is priced.

## 8. Prior art and reading limits

This is a finite model-specific certificate and experiment-language comparison,
not a discovery of nonlinear response or bilinear Hankel theory.

1. Petreczky, Wisniewski and Leth, *Moment matching for bilinear systems with nice
   selections*, arXiv:1605.04414v1; IFAC-PapersOnLine 49(18), 838-843 (2016),
   DOI 10.1016/j.ifacol.2016.10.270. PRIMARY_TEXT_READ for sections 2.1-2.2 and
   3: homogeneous bilinear form, word coefficients, reachability/observability
   and minimum realizations. https://arxiv.org/html/1605.04414v1
2. Lucarini, *Interpretable and Equation-Free Response Theory for Complex
   Systems*, arXiv:2502.07908. PRIMARY_TEXT_READ for the retrieved v1 HTML,
   section II.2, equations (13)-(19): Markov-chain second-order response,
   time-ordering and two-mode spectral weights. Publisher indexed text was also
   retrieved (DOI 10.1098/rsta.2025.0081), but direct DOI opening returned 403;
   journal equation numbering is not substituted for the arXiv numbering.
   https://arxiv.org/html/2502.07908v1
3. Mueller, Basu, Sollich and Krueger, *Coarse-grained second-order response
   theory*, Physical Review Research 2, 043123 (2020),
   DOI 10.1103/PhysRevResearch.2.043123. ABSTRACT_ONLY plus selected publisher
   excerpts: an equilibrium second-order coarse-grained framework and its
   non-Markovian issues. Its equilibrium hypotheses are not silently imported
   to P398. https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.2.043123

These are bounded primary readings, not a literature-completeness/novelty
search. The older inference "no named Markov selection rule was found, hence
our statement is new" is unsupported: absence of a retrieved name is not an
originality certificate. The concrete P398 ranks and scalar observer above
remain the result to evaluate on their own merits.

## 9. Decision and reproduction

    python scripts/p398_double_pulse_exact.py --out /tmp/p398-pulses.json
    python scripts/verify_p398_double_pulse_certificate.py
    python -m unittest discover -s tests -p 'test_p398_double_pulse_exact.py'
    python scripts/p398_double_pulse_numerical_control.py

The first three commands need only the standard library. The last needs mpmath.
Output paths use exclusive creation; existing artifacts are not overwritten.

This completes a bounded second-order *accessibility* question in #594/#598/
#610, rather than commissioning more widths or repeating the old quadratic
amplitude ladder. It does not settle P398's approximate noisy mode-resolution
cost, an all-width growth law, or any square-site original-U map in #275.
A further large run requires an explicit alternative and signal/noise objective;
there is no such justified purchase in this delivery.
