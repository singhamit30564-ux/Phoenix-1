# Scientific audit of Phoenix

This document records the audit performed on the Phoenix program suite before
the scientific corrections in this branch. It is written against the PHYLIP 3.69
documentation (`phylipweb.github.io/phylip/doc/`), not against memory.

The audit question for every program was: **is the algorithm that runs the one
the PHYLIP name denotes?** A program name appearing in an interface is not
evidence of anything.

## Summary of findings

| | Count |
| --- | --- |
| Programs whose algorithm matched their PHYLIP counterpart | 25 |
| Programs whose algorithm did **not** match and were corrected in this branch | 4 |
| Programs whose algorithm partially matched and are now labelled | 6 |

Four mappings were wrong and are corrected here:

1. **DNAMOVE** was implemented as Camin–Sokal parsimony. It is *interactive DNA
   parsimony*.
2. **MOVE** was implemented as Camin–Sokal parsimony. It is *interactive mixed
   parsimony*.
3. **FACTOR** was implemented as principal coordinate analysis. It is
   *multistate-to-binary recoding*.
4. **CONTML** was implemented as least squares on Brownian-motion distances. It
   is *restricted maximum likelihood under Brownian motion*.

A fifth program was corrected for a subtler reason:

5. **DNAINVAR** counted site patterns supporting each split. PHYLIP implements
   Cavender's and Lake's *invariants*. A pattern count is not an invariant test,
   so the implementation was replaced with a genuine rank-condition
   computation.

## Per-program audit

Legend for **Status**: *validated* (checked numerically in
`validation/validate.js`), *documented differences*, *partial*, *exploratory*.
This table is reproduced from the live registry by
`validation/generate-compatibility.js`; see `docs/COMPATIBILITY.md` for the
generated version.

### Parsimony and compatibility

| Program | PHYLIP purpose | Implementation before this branch | Discrepancy | Correction |
| --- | --- | --- | --- | --- |
| PROTPARS | Protein parsimony | Fitch unordered parsimony | PHYLIP uses a constrained amino-acid step matrix; ours is unordered | Kept, labelled *documented differences* |
| DNAPARS | DNA parsimony | Fitch parsimony, stepwise + NNI | None material | None |
| DNAPENNY | DNA parsimony by branch and bound | Heuristic multi-restart NNI, described as "branch and bound" | **Misleading description**: no bounds are computed | Description corrected; labelled *partial* |
| DNAMOVE | **Interactive DNA parsimony** | **Camin–Sokal parsimony** | **Wrong algorithm** | Replaced with an interactive session: initial tree, per-site reconstructions, NNI rearrangement list with score deltas, undo |
| DNACOMP | DNA compatibility | Largest compatible site clique, tree from it | None material | None |
| PARS | Wagner parsimony for multistate characters | Sankoff with ordered costs | None material | None |
| MIX | Mixed Wagner / Dollo / Camin–Sokal parsimony | Weighted sum of the three criteria | None material | None |
| PENNY | Wagner parsimony by branch and bound | Heuristic multi-restart search, described as branch and bound | **Misleading description** | Description corrected; labelled *partial* |
| MOVE | **Interactive mixed parsimony** | **Camin–Sokal parsimony** | **Wrong algorithm** | Replaced with an interactive session using the mixed criterion with user weights |
| DOLLOP | Dollo parsimony | Dollo DP: one gain, unlimited losses | Multistate states were read only as 0 vs 1 | Non-zero states now read as "present", documented |
| DOLPENNY | Dollo by branch and bound | Heuristic multi-restart search | Description overstated the search | Labelled *partial* |
| DOLMOVE | Dollo + Camin–Sokal parsimony | Weighted sum of both criteria | None material | None |
| CLIQUE | Character compatibility | Estabrook compatibility + maximum clique | None material | None |

### Likelihood

| Program | PHYLIP purpose | Implementation before this branch | Discrepancy | Correction |
| --- | --- | --- | --- | --- |
| DNAML | DNA maximum likelihood | Felsenstein pruning, NNI, branch optimisation | Correct in structure; now numerically validated against brute-force enumeration | Validation added |
| DNAMLK | DNA ML with a clock | Constrained node-height optimisation | Optimiser differs from PHYLIP's; stated | Labelled *documented differences* |
| PROML | Protein ML | Poisson / equal-input rates | **Missing** empirical exchangeability matrices (JTT, Dayhoff, WAG, VT, PMB) | Documented; **not** claimed as supported |
| PROMLK | Protein ML with clock | As PROML with a clock | As above | Documented |
| RESTML | Restriction-site ML | Two-state likelihood | Correct in structure | Labelled *documented differences* |

### Distances

| Program | PHYLIP purpose | Implementation before this branch | Discrepancy | Correction |
| --- | --- | --- | --- | --- |
| DNADIST | DNA distances | p, JC, K2P closed forms; other models by numerical ML | Correct; JC and K2P validated against closed forms | Validation added |
| PROTDIST | Protein distances | p, Poisson, equal-input, ML | No empirical matrices | Documented |
| RESTDIST | Restriction-site distances | fraction, JC-style, ML | None material | None |
| GENDIST | Gene-frequency distances | Nei, Cavalli-Sforza chord, Reynolds | Nei's minimum distance not offered | Documented |
| TREEDIST | Tree distances | Robinson–Foulds + branch score | Quartet distance not offered | Documented |

### Distance trees

| Program | PHYLIP purpose | Implementation before this branch | Discrepancy | Correction |
| --- | --- | --- | --- | --- |
| NEIGHBOR | Neighbor joining / UPGMA | NJ and UPGMA | None material | Validation added (exact recovery on additive matrices) |
| FITCH | Fitch–Margoliash least squares | Weighted least squares + NNI | None material | Validation added |
| KITSCH | Fitch–Margoliash with clock | Ultrametric least squares | None material | Validation added |
| CONTML | **Restricted ML under Brownian motion** | **Least squares on −log correlation distances**, described as an approximation | **Wrong method**: no likelihood was computed | Replaced with restricted (REML) likelihood, root- and scale-invariance verified, density cross-checked independently |
| CONTRAST | Independent contrasts | Contrasts computed | PHYLIP also estimates variance components | Labelled *partial* |

### Resampling, diagnostics, drawing

| Program | PHYLIP purpose | Implementation before this branch | Discrepancy | Correction |
| --- | --- | --- | --- | --- |
| SEQBOOT | Bootstrap/jackknife resampling | Column resampling, seeded | Species resampling for frequency data not offered | Documented |
| CONSENSE | Consensus tree | Strict, majority-rule, custom threshold | None material | Validation added |
| RETREE | Tree editing | Reroot, resolve, rearrange, clean | Interactive branch-length editing missing | Labelled *partial* |
| DNAINVAR | **Cavender / Lake invariants** | **Per-site pattern counts supporting each split** | **Not an invariant method** | Replaced with Cavender's quadratic invariants as a rank condition on the reduced pattern matrix; simulation-validated |
| FACTOR | **Multistate → binary recoding** | **Principal coordinate analysis** | **Unrelated algorithm** | Replaced with character-state-tree recoding, reversible, with a round-trip check |
| DRAWGRAM | Tree drawing | Rectangular phylogram / cladogram | Many PHYLIP styles missing | Labelled *partial* |
| DRAWTREE | Unrooted tree drawing | Radial phylogram | Styles and formats missing | Labelled *partial* |

## Shared-code findings

- **Duplicated logic.** `characterStates` was defined twice in
  `phoenix-support.js`; the duplicate is removed.
- **Numerical stability.** The likelihood evaluator recomputed the matrix
  exponential once per site per edge. It now computes one transition matrix per
  edge per rate class and applies it to all sites: identical results, roughly an
  order of magnitude faster, which is what makes bootstrap over likelihood
  practical.
- **Reproducibility.** Every randomised procedure takes an explicit seed:
  `SEQBOOT`, `RETREE`, `CONTRAST` polytomy resolution and the interactive
  sessions. Validation checks that the same seed gives the same answer.
- **Input validation.** Duplicate names, unequal lengths, invalid characters,
  empty data and malformed Newick all raise typed errors with the offending
  value named. A column of missing data is verified not to change either the
  parsimony score or a pairwise distance.

## Known scientific limits after this branch

These are stated rather than hidden:

- Searches are heuristic. `DNAPENNY`, `PENNY` and `DOLPENNY` are exhaustive over
  addition orders, not true branch-and-bound.
- No empirical amino-acid exchangeability matrices are bundled. Protein
  likelihood and distances use Poisson or equal-input rates.
- `DNAINVAR` implements Cavender's quadratic invariants but not Lake's linear
  invariants, and performs no significance test.
- `CONTML` maximises the correct restricted likelihood by coordinate-wise
  golden-section search; Felsenstein's EM algorithm is not reproduced.
- `PROTPARS` uses the unordered Fitch criterion rather than PHYLIP's constrained
  step matrix.
- No alignment is performed; sequences must already be aligned.
