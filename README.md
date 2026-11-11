# Phoenix

A browser-based, PHYLIP-inspired phylogenetics workspace. Paste or upload aligned sequences, pick a program, and get a tree, a distance matrix or a report — with the result drawn as SVG and exportable as Newick. Everything runs locally in your browser; no sequences are sent to a server.

## Run

Open `index.html` in a modern browser, or serve the folder with `python3 -m http.server 8000` and visit `http://localhost:8000`. No build step or dependencies are required. The interface uses Google Fonts when available, with local font fallbacks.

## Two ways in

1. **Workspace** — the quick path: aligned DNA, a distance model, UPGMA, a tree.
2. **Program suite** — the full toolchest. Pick any of the thirty-five PHYLIP programs from the sidebar, supply data, set options, run.

## The program suite

Every program reads FASTA or PHYLIP input (numeric and tree programs take their own formats), validates it, and reports errors without discarding your data.

### Parsimony

| Program | What it does |
| --- | --- |
| `DNAPARS` | DNA parsimony, Fitch criterion, stepwise addition + NNI |
| `DNAPENNY` | DNA parsimony with an intensive branch-and-bound style search |
| `DNAMOVE` | DNA parsimony under the Camin–Sokal (irreversible) assumption |
| `DNACOMP` | DNA compatibility: tree from the largest compatible site set |
| `PROTPARS` | Protein parsimony, unordered amino-acid replacements |
| `PARS` | Wagner parsimony for ordered multistate characters |
| `PENNY` | Wagner parsimony with intensive search |
| `MOVE` | Camin–Sokal parsimony (state 0 ancestral, no reversals) |
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
| `CONTML` | Tree from continuous characters or allele frequencies (Brownian motion) |
| `CONTRAST` | Felsenstein's independent contrasts on a supplied tree |

### Resampling, diagnostics, drawing

| Program | What it does |
| --- | --- |
| `SEQBOOT` | Bootstrap or jackknife replication, optionally searched and summarised into a supported tree |
| `CONSENSE` | Strict, majority-rule or threshold consensus from many trees |
| `RETREE` | Reroot on an outgroup, resolve polytomies, rearrange branches |
| `DNAINVAR` | Per-site invariant counts for the three quartet splits |
| `FACTOR` | Principal coordinate analysis of a distance matrix |
| `DRAWGRAM` | Rectangular phylogram or cladogram |
| `DRAWTREE` | Radial unrooted phylogram |

## Input formats

- **Sequences** — FASTA, or sequential/interleaved PHYLIP (`<taxa> <sites>` on the first line).
- **Discrete characters** — same formats, states `0`–`9`.
- **Restriction sites** — same formats, `1` present and `0` absent.
- **Numeric data** — optional `taxa characters` header, then `name value value …` per line.
- **Trees** — Newick, one tree per line for multi-tree programs.

Missing data is handled per method: unobserved states are skipped when a distance is computed and treated as compatible with anything in parsimony.

## Fidelity and limitations

Phoenix is an independent reimplementation, not a port of the PHYLIP source. Read this before treating any number as authoritative.

- **Search is heuristic.** Parsimony and likelihood searches use stepwise addition, nearest-neighbour interchange and random restarts. They find good trees, not guaranteed optima. `DNAPENNY`, `PENNY`, `DOLPENNY` and `DNAMOVE` say "branch and bound" in PHYLIP's sense; here they run many more independent searches, which is exhaustive only for very small data sets.
- **Protein models are limited.** Empirical exchangeability matrices (JTT, Dayhoff, WAG, VT, PMB) are not bundled. Protein likelihood and distance use Poisson or equal-input rates, and the interface says so.
- **CONTML** uses Brownian-motion distances (`-log` correlation) followed by least-squares branch fitting, not a full gene-frequency likelihood.
- **CONTRAST** resolves polytomies arbitrarily with near-zero branches and warns you when it does.
- **DNAINVAR** counts site patterns supporting each quartet split. It applies no significance correction and works on exactly four taxa.
- **No alignment.** Sequences must already be aligned and of equal length. Phoenix does not align, translate or edit sequences.
- **Scale.** The interface is tuned for tens of taxa, not thousands. Likelihood on more than a few dozen taxa will feel slow.
- **Not for clinical use.** Results are for exploration and teaching.

## Tests

Run `node --test tests/*.test.js` with Node.js 18 or newer. The suite covers parsing, distance models, tree building, parsimony scoring, likelihood, consensus, resampling and a run of all thirty-five programs.

## Layout

| File | Contents |
| --- | --- |
| `analysis.js` | The original small UPGMA engine used by the quick workspace |
| `phoenix-core.js` | Parsing, character encoding, Newick, tree primitives |
| `phoenix-models.js` | Substitution models, transition probabilities, pairwise distances |
| `phoenix-infer.js` | UPGMA, neighbor joining, least squares, parsimony, likelihood |
| `phoenix-support.js` | Bootstrap, consensus, tree comparison, PCoA, compatibility, drawing |
| `phoenix-programs.js` | The catalogue: thirty-five programs mapped onto the engine |
| `programs.js` | Program browser UI |
| `app.js`, `index.html`, `styles.css` | The workspace interface |

Each engine file uses a small UMD wrapper, so the same code runs in the browser and under Node for tests.

## Credits

Inspired by [PHYLIP](https://phylipweb.github.io/phylip/), created by Joseph Felsenstein. Phoenix is an independent reimplementation of a subset of phylogenetic functionality; it does not bundle PHYLIP or claim bit-for-bit compatibility with its programs.
