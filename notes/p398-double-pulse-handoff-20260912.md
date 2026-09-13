# Existing-issue update text: completed P398 double-pulse analysis

Use in #594 (with references to #598/#600/#610), not as a duplicate new task.
This text has not been posted to GitHub by this session.

---

Completed a bounded exact analysis using the ORIGINAL four sources and three
primary readouts of P398, widths 4/5. The change is the allowed intervention
time program, not a new source/readout or a reinterpretation of the old rank-6
approximation.

With H=join@0-join@(w-2) (twice archived H_odd), baseline G commutes with R,
H is odd, and G+epsilon H is physical for |epsilon|<=1. First-order even-to-even
response vanishes, as already known. The new result is that the two-insertion
kernel S H exp(tau G) H F has exact minimal delay order 4/16, exhausting the
odd subspaces. Finite positive pulse durations preserve these ranks for all
sufficiently small durations by analyticity of a nonzero Hankel minor.

Exact raw baseline order is 10/26. The full noncommutative G/H word response
has minimum homogeneous bilinear order 14/42; original source-minus-uniform
contrasts give 13/41 because the constant mode is invisible. Every rank has an
explicit nonzero integer Hankel minor and an independent mathematical upper
bound; all ten determinant witnesses were rechecked with integer Bareiss
elimination, including reconstruction of their G/H words.

At width 4, the old wrapped-pair source and wrap readout already give

    Khat(z)=2(z^2+11z+27)/[(z^2+10z+23)(z^2+11z+26)],
    K(0)=0, K'(0)=2, Hankel determinant=-16.

Thus zero delay can miss a fully observable four-mode delay response. Exact
opposite-sign examples in the old dictionary also show that parity fixes no
universal sign of the quadratic term.

Files: `notes/p398-double-pulse-observability-20260912.md`, the exact generator/
certificate reader and independent verifier, a finite-window numerical control,
and additive JSON results. Seven local mathematical tests pass. No new MC,
full-repository CI, large widths, GPU, physical p_c or continuum claim. Generic
response and bilinear realization formulas have primary prior art; the result
is the finite P398 observability certificate, not a claim that those tools are new.

This completes the small second-order accessibility question. No default width
extension is requested. A noise-limited identification task, if later justified,
needs an explicit alternative and error budget; exact full rank alone is not a
sample-size plan. #275's original-U candidate-map requirement remains separate.
