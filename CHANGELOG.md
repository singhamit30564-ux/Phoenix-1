# Changelog

## Scientific audit and validation (branch `feature/scientific-audit-and-validation`)

### Corrected algorithm mappings

- **DNAMOVE** is now interactive DNA parsimony, matching PHYLIP: it starts from a
  tree, reports the most-parsimonious reconstruction of each site, and lets you
  rearrange the topology by hand while the step count updates. It was previously
  (incorrectly) Camin–Sokal parsimony.
- **MOVE** is now interactive mixed parsimony using the Wagner / Dollo /
  Camin–Sokal weights of MIX. It was previously (incorrectly) Camin–Sokal
  parsimony.
- **FACTOR** now recodes multistate characters into binary characters using a
  character-state tree, with a reversible-encoding check and optional explicit
  state trees. It was previously (incorrectly) principal coordinate analysis.
- **CONTML** now maximises the Brownian-motion restricted (REML) likelihood,
  with root- and scale-invariance verified. It was previously least squares on
  −log correlation distances.
- **DNAINVAR** now computes Cavender's quadratic invariants as a rank condition
  on the reduced pattern matrix. It previously counted supporting site patterns,
  which is not an invariant method.

### Additions

- `phoenix-linalg.js` — Cholesky factorisation, inverse, log-determinant.
- `phoenix-continuous.js` — Brownian-motion restricted likelihood and tree
  fitting.
- `phoenix-interactive.js` — interactive parsimony sessions with undo, greedy
  completion and per-site reconstructions.
- `validation/validate.js` — 36 checks against independent references, writing
  `validation/report.json`.
- `validation/generate-compatibility.js` — generates `docs/COMPATIBILITY.md`
  from the live registry.
- `docs/AUDIT.md`, `docs/METHODS.md`, `docs/COMPATIBILITY.md`.
- Every program now declares a scientific status and a reference note, shown in
  the interface as a badge.

### Fixes

- Removed a duplicated `characterStates` definition.
- Likelihood evaluation now computes one transition matrix per edge per rate
  class instead of one per site, giving the same results roughly ten times
  faster.
- Dollo parsimony now treats any non-zero state as "present" for multistate data
  instead of reading only state 1.
- `DNAPENNY`, `PENNY`, `DOLPENNY` no longer describe themselves as
  branch-and-bound, which they are not.
- CONTML reports branch lengths in units of the Brownian-motion rate and states
  that only relative lengths are identifiable.
- The interface now distinguishes interactive programs and shows the
  validation status of each program.

### Known limitations (unchanged, now documented)

- Heuristic searches; no branch-and-bound.
- No empirical amino-acid exchangeability matrices (JTT, Dayhoff, WAG).
- DNAINVAR implements Cavender's invariants only, with no significance test.
- CONTML uses coordinate-wise optimisation rather than Felsenstein's EM.
- PROTPARS uses the unordered Fitch criterion.
- No alignment is performed.

## Program suite (branch `feature/phylip-program-suite`)

- Added all 35 PHYLIP programs, the engine modules `phoenix-core.js`,
  `phoenix-models.js`, `phoenix-infer.js`, `phoenix-support.js`,
  `phoenix-programs.js`, the program browser UI, and an initial test suite.
