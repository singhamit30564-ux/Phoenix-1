# PHYLIP compatibility matrix

This file is generated from the program registry by
`node validation/generate-compatibility.js`. Do not edit it by hand; edit
`phoenix-programs.js` and regenerate.

Phoenix is an independent reimplementation. It is **not** a port of the PHYLIP
source, and no program here claims bit-for-bit agreement with PHYLIP output.
Status labels mean:

- **validated** — Fully implemented, with numerical checks against an independent reference.
- **documented differences** — Correct method, but a feature or optimisation of the original differs and is documented.
- **partial** — Part of the original behaviour is missing.
- **exploratory** — A reasoned approach whose numbers should not be treated as reference values.

## Summary

| Status | Programs |
| --- | --- |
| validated | 15 |
| documented differences | 12 |
| partial | 8 |
| exploratory | 0 |
| **total** | **35** |

Only a program whose method is checked numerically by
`validation/validate.js` is labelled **validated**.

### Parsimony and compatibility

| Program | Status | What is implemented | Difference from PHYLIP |
| --- | --- | --- | --- |
| `PROTPARS` | documented differences | Protein parsimony under the Fitch criterion, treating amino-acid replacements as unordered and equally costly. | Unordered Fitch criterion; PHYLIP PROTPARS uses a constrained amino-acid step matrix |
| `DNAPARS` | validated | Finds the most parsimonious tree for DNA data using the Fitch (unordered) criterion, with stepwise addition and nearest-neighbour interchange. | Fitch (1971) parsimony; step counts checked against hand-computed trees |
| `DNAPENNY` | partial | DNA parsimony searched exhaustively over addition orders and nearest-neighbour rearrangements, keeping the shortest tree found. This is exact only for small data sets; PHYLIP uses branch-and-bound with bounds and can go further. | Exhaustive over addition orders and NNI; PHYLIP uses branch-and-bound with bounds |
| `DNAMOVE` | documented differences | Interactive DNA parsimony, after MacClade. Start from a tree, read the most-parsimonious reconstruction of every site, then rearrange the tree by hand and watch the step count respond. | Interactive parsimony as in PHYLIP, presented through a web interface instead of a terminal |
| `DNACOMP` | validated | Finds the largest set of mutually compatible sites (the clique of characters whose state distributions admit a single tree) and builds the tree implied by them. | Character compatibility; clique reported and checked |
| `PARS` | validated | Wagner parsimony for multistate discrete characters: states 0–9 are ordered, and a change of k states costs k steps. | Wagner parsimony with ordered states |
| `MIX` | validated | Weighs Wagner, Dollo and Camin–Sokal criteria against each other so that no single assumption about character evolution dominates. | Weighted Wagner / Dollo / Camin-Sokal sum, as in PHYLIP MIX |
| `PENNY` | partial | Wagner parsimony with an intensive search over addition orders and nearest-neighbour interchange. | Intensive heuristic search; PHYLIP uses branch-and-bound |
| `MOVE` | documented differences | Interactive parsimony for multistate characters using the mixed Wagner / Dollo / Camin–Sokal criterion of MIX. Rearrange the tree by hand and watch the weighted step count change. | Interactive mixed parsimony, web interface rather than terminal |
| `DOLLOP` | validated | Dollo parsimony: a derived state may be gained only once, but may be lost repeatedly. Suits presence/absence data such as restriction sites. | Dollo criterion: one gain, unlimited losses |
| `DOLPENNY` | partial | Dollo parsimony with intensive search over addition orders. | Intensive heuristic search; PHYLIP uses branch-and-bound |
| `DOLMOVE` | documented differences | Combines the Dollo criterion with the Camin–Sokal assumption for the two-state case, allowing one gain per character while barring reversals. | Weighted Dollo and Camin-Sokal combination |
| `CLIQUE` | validated | Reports the largest clique of pairwise compatible characters. Two characters are compatible when the states of one are nested or disjoint within the other. | Estabrook compatibility and maximum clique |

### Likelihood

| Program | Status | What is implemented | Difference from PHYLIP |
| --- | --- | --- | --- |
| `DNAML` | validated | Maximum-likelihood tree for DNA under a reversible model with optional discrete-gamma rate variation. Branch lengths are optimised and the topology refined by nearest-neighbour interchange. | Felsenstein (1981) pruning likelihood; checked by brute-force enumeration on four taxa |
| `DNAMLK` | documented differences | Maximum likelihood for DNA with a molecular clock enforced, so all tips are equidistant from the root. | Clock enforced by constrained node-height optimisation, not by PHYLIP's own search |
| `PROML` | documented differences | Maximum-likelihood tree for amino-acid sequences with Poisson or equal-input rates and optional gamma rate variation. | Poisson / equal-input rates; empirical exchangeability matrices are not bundled |
| `PROMLK` | documented differences | Protein maximum likelihood with a strict molecular clock. | As PROML, with a clock |
| `RESTML` | documented differences | Maximum likelihood for presence/absence restriction-site data using a two-state model. | Two-state likelihood; rate variation only through gamma categories |

### Distances

| Program | Status | What is implemented | Difference from PHYLIP |
| --- | --- | --- | --- |
| `DNADIST` | validated | Pairwise evolutionary distances between DNA sequences, using a choice of substitution models and analytic or maximum-likelihood corrections. | JC69, K2P and ML distances checked against closed forms |
| `PROTDIST` | documented differences | Pairwise distances between amino-acid sequences: observed proportion of differences, Poisson correction, equal-input correction or maximum likelihood. | No empirical exchangeability matrices; Poisson and equal-input only |
| `RESTDIST` | validated | Pairwise distances from restriction-site presence/absence data: the simple fraction of differing sites or a maximum-likelihood estimate under a two-state model. | Two-state corrections checked against closed forms |
| `GENDIST` | documented differences | Genetic distances between populations from allele frequencies: Nei’s standard distance, Cavalli-Sforza chord distance or Reynolds’ coancestry distance. | Nei, Cavalli-Sforza chord and Reynolds distances; Nei's minimum distance is not offered |
| `TREEDIST` | documented differences | Compares trees: the symmetric difference (Robinson–Foulds) counts splits present in one tree but not the other, and the branch score compares branch lengths. | Robinson-Foulds and branch-score distances; quartet distances are not offered |

### Distance trees

| Program | Status | What is implemented | Difference from PHYLIP |
| --- | --- | --- | --- |
| `NEIGHBOR` | validated | Builds a tree from a distance matrix using neighbour joining (Saitou–Nei) or UPGMA. Both run on the distances implied by the model you choose. | Saitou-Nei neighbor joining; exact recovery on additive matrices |
| `FITCH` | validated | Least-squares tree without a clock: branch lengths are chosen to minimise the weighted squared difference between observed and tree distances, and the topology is refined by NNI. | Fitch-Margoliash weighted least squares; recovers additive distances exactly |
| `KITSCH` | validated | Least-squares tree constrained to a molecular clock, so every tip is the same distance from the root. | Ultrametric least squares; ultrametricity verified |
| `CONTML` | documented differences | Restricted maximum likelihood under the Brownian-motion model of Cavalli-Sforza and Edwards, for continuous characters or allele frequencies. Branch lengths are fitted by maximising the restricted likelihood and the topology is refined by NNI. | Brownian-motion restricted ML with coordinate-wise optimisation; the EM algorithm is not reproduced |
| `CONTRAST` | partial | Computes Felsenstein’s independent contrasts for continuous characters on the supplied tree: each contrast is the difference between sister lineages, scaled by the variance expected under Brownian motion. | Independent contrasts computed; PHYLIP CONTRAST also estimates variance components |

### Resampling and diagnostics

| Program | Status | What is implemented | Difference from PHYLIP |
| --- | --- | --- | --- |
| `SEQBOOT` | validated | Generates bootstrap replicate data sets by resampling alignment columns with replacement, or summarises them straight into a supported tree with a consensus. | Column resampling with a reproducible seed; proportions checked on a known tree |
| `CONSENSE` | validated | Builds a consensus tree from a collection of trees: strict, majority rule, or a threshold you choose. | Strict and majority-rule consensus with support percentages |
| `RETREE` | partial | Reads a Newick tree and rewrites it: reroot on an outgroup, resolve polytomies, or apply a random rearrangement. | Reroot, resolve, rearrange, clean; interactive branch-length editing is not offered |
| `DNAINVAR` | partial | Computes Cavender's invariants for four taxa. Nucleotides are reduced to purines and pyrimidines, the 4×4 pattern matrix is formed for each of the three topologies, and the invariant residual measures how far each matrix is from the rank-two structure that topology requires. | Cavender quadratic invariants as a rank test; Lake linear invariants and significance tests are not implemented |
| `FACTOR` | validated | Recodes multistate characters into binary characters using a character-state tree, so the result can be analysed by the two-state programs. Reports the mapping and verifies that the recoding is reversible. | Character-state-tree recoding with a reversibility check |

### Drawing

| Program | Status | What is implemented | Difference from PHYLIP |
| --- | --- | --- | --- |
| `DRAWGRAM` | partial | Renders a rooted tree as a rectangular phylogram or cladogram, drawn directly in the browser. | Rectangular and cladogram styles only; PHYLIP offers many typographic styles |
| `DRAWTREE` | partial | Renders an unrooted tree as a radial phylogram, the classic DRAWTREE view. | Radial phylogram only; PHYLIP offers several styles and formats |
