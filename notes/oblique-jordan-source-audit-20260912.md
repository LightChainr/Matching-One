# Source and action audit for the oblique/Jordan delivery

2026-09-12. This is a bounded primary-source recheck, not a global novelty survey.

## Repository reads

- Latest PR search: #732/#731/#730/#729/#728/#725/#724/#723/#718. Full PR
  bodies were read, not taken as endorsements of each other's conclusions.
- #724 at ef6bdcf5c8f77fa697f0273731bb681334a1967e:
  notes/lit-thermal-jet-jordan-20260912.md read directly through GitHub.
- Actual #724 and #728 discussions read before adding corrections.
- #718's aspect-uniform note and rectangular reference script read from the
  conversation's supplied immutable archive. New scripts have NO runtime
  dependency on that archive or on #708/#710/#716/#718.
- #732's peer check supports the stated rectangular result. It does not
  certify the new arbitrary-period theorem, which needs its own review.
- #728 keeps a false distinction between Mertens--Ziff's scalar matching
  function and digital-Alexander M. The previous site-source handoff already
  corrected that distinction; this round posts it where it was propagated,
  rather than counting it as another new theorem.

## External primary text read this round

1. Duminil-Copin--Tassion, arXiv:1502.03050v3,
   https://arxiv.org/html/1502.03050v3
   Thm 1.1(3), section 1.2 site adaptation. PRIMARY_TEXT_READ for these sections.
   The printed theorem is bond language, the site paragraph is an adaptation.
   No new claim that the paper prints a separate numbered site theorem.
2. Grimmett--Li, arXiv:2205.02734v3,
   https://arxiv.org/html/2205.02734v3
   Introduction, especially Eq (1.3), matching-pair relation and amenable scope.
   PRIMARY_TEXT_READ for these sections. Companion proof 2203.00981 was not
   separately read this round and is not given a new primary-read tag.
3. Bamieh, arXiv:2002.05001v2,
   https://arxiv.org/html/2002.05001v2
   Assumptions and section 2: analytic eigenvectors/eigenvalues are explicitly
   assumed. PRIMARY_TEXT_READ. Kato's body is NOT read or quoted this round.
4. Vasseur--Jacobsen--Saleur, arXiv:1206.2312v2,
   https://arxiv.org/html/1206.2312v2
   Eq (11) and conclusion: specific field mixing vs Boltzmann-weight derivative
   logs. PRIMARY_TEXT_READ for these passages. No general novelty of the finite
   matrix counterexamples or visibility expansion is claimed.
5. Mertens--Ziff, arXiv:1603.07289v2,
   https://arxiv.org/html/1603.07289v2
   Introduction after (4), (20)--(21) and following text. PRIMARY_TEXT_READ for
   these passages. They explicitly connect the matching RHS and critical-polynomial
   criterion. Same scalar object does not identify a particular transfer matrix.
6. Qian--Chu--Tan, A Systematic Analysis on Analyticity of Semisimple Eigenvalues
   of Matrix-Valued Functions, DOI 10.1137/15M1053050.
   Publisher abstract read; ABSTRACT_ONLY. Not used as the proof of a counterexample.

No PDF figures or tables were needed; the above reads used primary HTML.

## Review decisions

- Keep the existing site-sharpness/matching assumptions and #718 rectangular
  conclusion. Extend them using an explicit new oblique-band proof, not an
  unproved change of physical nearest-neighbour geometry.
- Correct #724/#731's linear-split implication and analytic-branch assumption.
  Keep the narrow warning that a thermal polynomial length factor is not enough
  to identify a physical Jordan block.
- Correct #728's scalar-object distinction; keep the absence of our specific
  small block representation in those cited equations.
- Bounded retrieval that does not find a formula should be stated as such.
  Claims in #723/#730 of absence from all print / surviving novelty should not
  be upgraded to priority certification by this delivery.

## Remote actions actually performed

Top-level comments were posted successfully:

- #724, comment 5646908635: exact Jordan/analyticity counterexamples and visibility.
- #731, comment 5646909623: peer-check amendment, referring to the #724 correction.
- #728, comment 5646913116: same scalar observable vs different matrix realization.
- #718, comment 5646924543: completed arbitrary-period proof and marked-source scope.

No branch writes, merge, issue creation/closure, STATUS edit, workflow trigger,
production request or credential/hardware operation occurred. The executable
files are handed to the owner as an additive patch; comments are not substitutes
for committing or independently reviewing the proof.
