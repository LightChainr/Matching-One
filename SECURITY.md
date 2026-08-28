# Security Policy

Matching One is research software, not a production service. Security reports are nevertheless important when they concern repository integrity, dependency compromise, arbitrary code execution, malicious input handling, CI credentials, or falsification of research provenance.

## Supported code

Security fixes target the current `main` branch and active integration pull requests. Historical result branches and immutable archives may receive a warning or correction record rather than an in-place patch.

## Reporting

Do not open a public issue containing an exploit, credential, private token, or weaponized proof of concept.

Use GitHub's private security-advisory or private vulnerability-reporting interface for this repository when available. If that interface is unavailable, contact the maintainer through the GitHub profile and provide only enough public information to establish a private channel.

Include:

- affected commit and file;
- impact and threat model;
- reproduction steps or minimal proof of concept;
- whether secrets, CI permissions, released artifacts, or result provenance are affected;
- suggested mitigation if known.

## Research-integrity incidents

A defect that changes numerical results, covariance, RNG domains, geometry conventions, or provenance should normally be reported as a research-correction issue rather than a security vulnerability, unless exploitation or intentional tampering is involved.

Committed results are not silently rewritten. Corrections follow `REPRODUCIBILITY.md` and update `docs/STATUS.md` when claim strength changes.
