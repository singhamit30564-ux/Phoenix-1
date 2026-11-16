# Changelog

## Review remediation (PR #3 follow-up)

Fixes four defects found in the corrections themselves.

### MOVE scored with Fitch instead of Wagner

- The weighted branch in `phoenix-interactive.js` was guarded by
  `if (weights.dollo || weights.caminSokal)`, so the default configuration
  (Wagner 1, Dollo 0, Camin–Sokal 0) fell through to `fitchScore`. On a data set
  where the two criteria differ (`0,0,2,2` costs 2 under ordered Wagner and 1
  under unordered Fitch) the session now reports 2.
- `resolveCriterion()` replaces the fall-through: DNA and protein data use
  Fitch; discrete data resolve to exactly one named criterion or to `mixed`.
  All-zero or negative weights raise an error rather than degrading silently.
- An explicit `{ dollo: 1 }` is no longer widened to "Dollo plus Wagner"; only a
  completely absent weight specification defaults to Wagner.

### Reconstruction ignored the scored criterion

- `siteReconstruction()` used Fitch set intersection for every criterion. It now
  runs the same dynamic program over (gains used, state) that produces the
  score, so per-character reconstruction costs sum exactly to the reported
  score, in Wagner, Dollo, Camin–Sokal and mixed modes.

### FACTOR constant-character option did nothing

- `factorRecode()` dropped constant characters before `deleteConstant` was
  consulted, so the option had no effect. The option is removed.
- Behaviour now follows PHYLIP: with any state tree declared, characters without
  one are omitted; a character whose states are all identical has a one-node
  tree, no edges, and contributes no binary characters. Both cases are reported
  separately.

### Character-state trees were not validated

- `validateStateTree()` rejects self-loops, repeated edges, cycles,
  disconnection, invalid identifiers and observed states missing from the tree,
  naming the character in every message.

### Documentation

- `docs/AUDIT.md` records the second remediation round and reconciles the status
  counts with the program registry (15 validated, 12 documented differences,
  8 partial, 0 exploratory).
- `docs/METHODS.md` documents the criterion dispatch, the reconstruction dynamic
  program and the FACTOR character-handling rules.
- `validation/PHYLIP-DIFFERENTIAL.md` states explicitly that no PHYLIP
  differential comparison has been performed, and records the protocol for one.
- `validation/validate.js` grows from 36 to 58 checks, including a check that the
  compatibility matrix counts match the registry.


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
