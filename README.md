# Phoenix

A browser-based, PHYLIP-inspired phylogenetics workspace. Paste or upload aligned sequences, pick a program, and get a tree, a distance matrix or a report — with the result drawn as SVG and exportable as Newick. Everything runs locally in your browser; no sequences are sent to a server.

## Run

Open `index.html` in a modern browser, or serve the folder with `python3 -m http.server 8000` and visit `http://localhost:8000`. No build step or dependencies are required. The interface uses Google Fonts when available, with local font fallbacks.

## Two ways in

1. **Workspace** — the quick path: aligned DNA, a distance model, UPGMA, a tree.
2. **Program suite** — the full toolchest. Pick any of the thirty-five PHYLIP programs from the sidebar, supply data, set options, run. Each program shows a badge saying how far its implementation has been validated.

## The program suite

Every program reads FASTA or PHYLIP input (numeric and tree programs take their own formats), validates it, and reports errors without discarding your data.

### Parsimony

| Program | What it does |
| --- | --- |
| `DNAPARS` | DNA parsimony, Fitch criterion, stepwise addition + NNI |
| `DNAPENNY` | DNA parsimony with an intensive branch-and-bound style search |
| `DNAMOVE` | Interactive DNA parsimony: rearrange the tree by hand, watch the steps change |
| `DNACOMP` | DNA compatibility: tree from the largest compatible site set |
| `PROTPARS` | Protein parsimony, unordered amino-acid replacements |
| `PARS` | Wagner parsimony for ordered multistate characters |
| `PENNY` | Wagner parsimony with intensive search |
| `MOVE` | Interactive mixed parsimony (Wagner / Dollo / Camin–Sokal weights) |
| `MIX` | Weighted combination of Wagner, Dollo and Camin–Sokal |
| `DOLLOP` | Dollo parsimony: one gain per character, unlimited losses |
| `DOLPENNY` | Dollo parsimony with intensive search |
| `DOLMOVE` | Dollo combined with Camin–Sokal |
| `CLIQUE` | Largest clique of mutually compatible characters |

### Likelihood

| Program | What it does |
| --- | --- |
| `DNAML` | DNA maximum likelihood (JC, F81, K2P, HKY, TN93, GTR; optional Γ rates) |
| `DNAMLK` | DNA maximum likelihood under a strict molecular clock |
| `PROML` | Protein maximum likelihood (Poisson / equal-input; optional Γ) |
| `PROMLK` | Protein maximum likelihood under a strict clock |
| `RESTML` | Restriction-site maximum likelihood (two-state model) |

### Distances

| Program | What it does |
| --- | --- |
| `DNADIST` | Pairwise DNA distances (p, JC, K2P, F81, HKY, TN93, GTR or ML) |
| `PROTDIST` | Pairwise protein distances (p, Poisson, equal-input or ML) |
| `RESTDIST` | Restriction-site distances (fraction, JC-style or ML) |
| `GENDIST` | Gene-frequency distances (Nei, Cavalli-Sforza chord, Reynolds) |
| `TREEDIST` | Robinson–Foulds and branch-score distances between trees |

### Distance trees

| Program | What it does |
| --- | --- |
| `NEIGHBOR` | Neighbor joining or UPGMA from a distance matrix |
| `FITCH` | Fitch–Margoliash least squares, no clock |
| `KITSCH` | Fitch–Margoliash least squares with a clock |
| `CONTML` | Brownian-motion restricted maximum likelihood for continuous characters and gene frequencies |
| `CONTRAST` | Felsenstein's independent contrasts on a supplied tree |

### Resampling, diagnostics, drawing

| Program | What it does |
| --- | --- |
| `SEQBOOT` | Bootstrap or jackknife replication, optionally searched and summarised into a supported tree |
| `CONSENSE` | Strict, majority-rule or threshold consensus from many trees |
| `RETREE` | Reroot on an outgroup, resolve polytomies, rearrange branches |
| `DNAINVAR` | Cavender's invariants for the three quartet splits |
| `FACTOR` | Multistate to binary recoding through a character-state tree |
| `DRAWGRAM` | Rectangular phylogram or cladogram |
| `DRAWTREE` | Radial unrooted phylogram |

## Input formats

- **Sequences** — FASTA, or sequential/interleaved PHYLIP (`<taxa> <sites>` on the first line).
- **Discrete characters** — same formats, states `0`–`9`.
- **Restriction sites** — same formats, `1` present and `0` absent.
- **Numeric data** — optional `taxa characters` header, then `name value value …` per line.
- **Trees** — Newick, one tree per line for multi-tree programs.

Missing data is handled per method: unobserved states are skipped when a distance is computed and treated as compatible with anything in parsimony.

## Scientific status

Phoenix is an independent reimplementation, not a port of the PHYLIP source, and it does not claim bit-for-bit agreement with PHYLIP. Every program carries a status:

| Status | Meaning |
| --- | --- |
| **validated** | The method is implemented and its numbers are checked against an independent reference by `validation/validate.js`. |
| **documented differences** | The method is the right one, but a feature or optimiser of the original differs, and the difference is written down. |
| **partial** | Part of the original behaviour is missing. |

The full matrix is in [`docs/COMPATIBILITY.md`](docs/COMPATIBILITY.md), generated from the program registry so it cannot drift from the code. The audit that led to the current state, including the four programs whose algorithms were wrong, is in [`docs/AUDIT.md`](docs/AUDIT.md). The mathematics is in [`docs/METHODS.md`](docs/METHODS.md).

## Validation

```bash
npm test            # 47 unit and integration tests
npm run validate    # 36 numerical checks, writes validation/report.json
npm run check       # both
```

`validation/validate.js` compares Phoenix output with values obtained independently — closed-form distances, exhaustive enumeration of internal states for the likelihood, hand-computed parsimony scores, additive matrices for neighbor joining, simulations with a known generating topology for the invariant test, and an independent Gauss–Jordan evaluation for the Brownian-motion likelihood. Each check records its expected value, the actual value, the tolerance and the reference it used, in machine-readable form.

These checks establish implementation correctness against mathematical references. They are not a validation against PHYLIP's own output, and they are not a substitute for comparing against published analyses.

## Fidelity and limitations

Read this before treating any number as authoritative.

- **Search is heuristic.** Parsimony and likelihood searches use stepwise addition, nearest-neighbour interchange and random restarts. They find good trees, not guaranteed optima. `DNAPENNY`, `PENNY` and `DOLPENNY` search exhaustively over addition orders but do not implement branch-and-bound with bounds, so they are exact only for small data sets.
- **Protein models are limited.** Empirical exchangeability matrices (JTT, Dayhoff, WAG, VT, PMB) are not bundled. Protein likelihood and distance use Poisson or equal-input rates, and the interface says so.
- **CONTML** maximises the correct Brownian-motion restricted likelihood, but only relative branch lengths are identifiable: the rate and the branch-length scale cannot both be estimated. Felsenstein's EM algorithm is not reproduced; Phoenix uses coordinate-wise golden-section search on the same criterion.
- **CONTRAST** resolves polytomies arbitrarily with near-zero branches and warns you when it does.
- **DNAINVAR** computes Cavender's quadratic invariants (a rank condition on the reduced pattern matrix). Lake's linear invariants and any significance test are not implemented; exactly four taxa are required.
- **FACTOR** recodes through a character-state tree derived from the observed states in numeric order, unless you declare one. PHYLIP's FACTOR takes the state trees from the input file.
- **No alignment.** Sequences must already be aligned and of equal length. Phoenix does not align, translate or edit sequences.
- **Scale.** The interface is tuned for tens of taxa, not thousands. Likelihood on more than a few dozen taxa will feel slow.
- **Not for clinical use.** Results are for exploration and teaching.

## Tests

Run `npm test` with Node.js 18 or newer. The suite covers parsing, distance models, tree building, parsimony scoring, likelihood, consensus, resampling, recoding, invariants and Brownian-motion likelihood, and runs all thirty-five programs.

## Layout

| File | Contents |
| --- | --- |
| `analysis.js` | The original small UPGMA engine used by the quick workspace |
| `phoenix-core.js` | Parsing, character encoding, Newick, tree primitives |
| `phoenix-models.js` | Substitution models, transition probabilities, pairwise distances |
| `phoenix-linalg.js` | Cholesky factorisation, inverse, log-determinant |
| `phoenix-infer.js` | UPGMA, neighbor joining, least squares, parsimony, likelihood |
| `phoenix-continuous.js` | Brownian-motion restricted likelihood |
| `phoenix-interactive.js` | Interactive parsimony sessions |
| `phoenix-support.js` | Bootstrap, consensus, tree comparison, PCoA, compatibility, drawing |
| `phoenix-programs.js` | The catalogue: thirty-five programs mapped onto the engine |
| `programs.js` | Program browser UI |
| `app.js`, `index.html`, `styles.css` | The workspace interface |
| `validation/` | Numerical validation suite and its machine-readable report |
| `docs/` | Audit, methods and compatibility matrix |

Each engine file uses a small UMD wrapper, so the same code runs in the browser and under Node for tests.

## Credits

Inspired by [PHYLIP](https://phylipweb.github.io/phylip/), created by Joseph Felsenstein. Phoenix is an independent reimplementation of a subset of phylogenetic functionality; it does not bundle PHYLIP or claim bit-for-bit compatibility with its programs.
