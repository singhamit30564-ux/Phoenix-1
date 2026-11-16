# Phoenix

**A browser-based reimplementation of the PHYLIP phylogenetics toolchest — with an explicit scientific status on every one of its 35 programs.**

Phoenix turns aligned sequence and character data into trees, distance matrices and reports entirely inside the browser. No install, no server, no data leaving the machine. Everything you see is computed by the code in this repository; nothing is simulated for display.

[Run it locally](#run-it) · [Program catalogue](#program-catalogue) · [Validation](#validation) · [Guided demos](#guided-demos) · [Known limitations](#known-limitations)

---

## Why this exists

PHYLIP has been the standard teaching and research package for phylogenetics since 1980. It is powerful, but it is a command-line package with a menu-driven interface that most students today have never used, and it is hard to demonstrate live.

Phoenix keeps the science and replaces the delivery: the same established algorithms, in a workspace you can open in a tab, with the input validated, the result explained, and every output exportable. It is an independent reimplementation — **not** a port of the PHYLIP source — and it says so everywhere it matters.

## Capabilities

- **35 programs** across parsimony, likelihood, distances, distance trees, resampling, diagnostics and tree drawing.
- **Explicit scientific status** on every program: validated, documented difference, or partial. Nothing is labelled "validated" without a numerical check behind it.
- **Real numerical work**: Felsenstein pruning with discrete-gamma rates, Fitch/Sankoff/Dollo/Camin–Sokal parsimony, Brownian-motion restricted likelihood, Cavender invariants, character-state-tree recoding.
- **Interactive parsimony** (DNAMOVE/MOVE): rearrange a tree by hand and watch the step count respond.
- **Reproducibility**: every randomised procedure takes a seed, and each result exports a JSON record of the program, parameters, input summary and outcome.
- **No dependencies.** One HTML file, plain scripts, no build step, no framework.

## Run it

The application is a static site. A local server is recommended because the Web Worker that keeps long analyses off the main thread is blocked on `file://`.

```bash
git clone https://github.com/singhamit30564-ux/Phoenix-1.git
cd Phoenix-1
python3 -m http.server 8000
# open http://localhost:8000
```

Node.js 18+ is only needed for the test and validation suites:

```bash
npm test              # 88 unit and integration tests
npm run validate      # 62 numerical checks against independent references
npm run compat        # regenerate the compatibility matrix from the registry
npm run check         # tests + validation
```

Optional, for browser checks and screenshots (requires Playwright):

```bash
npm i -g playwright && npx playwright install chromium
npm run browser-check   # 12 checks: all 35 programs, demos, exports, responsive layouts
npm run screenshots     # writes docs/screenshots/*.png
npm run benchmark       # writes validation/benchmark.json
```

## Guided demos

Three one-click workflows, each running the real programs and showing every intermediate step.

| Demo | What it does | Programs exercised |
| --- | --- | --- |
| **A — DNA distances and a tree** | Jukes–Cantor distances on five primate sequences, then neighbour joining | `DNADIST` → `NEIGHBOR` |
| **B — Parsimony and ancestral states** | Wagner parsimony on ordered characters, with the most-parsimonious reconstruction of each character | `PARS` → `MOVE` |
| **C — Bootstrap support** | 50 bootstrap replicates searched with parsimony, summarised by majority-rule consensus | `SEQBOOT` → `CONSENSE` |

Demo C uses sequences simulated under JC69 on a **known** tree, so the recovered clades (all at 100% support) can be checked against the truth.

## Program catalogue

Status is generated from the program registry, so this table cannot drift from the code. See [`docs/COMPATIBILITY.md`](docs/COMPATIBILITY.md) for the full matrix with per-program differences.

### Parsimony and compatibility

| Program | Status | What it does |
| --- | --- | --- |
| `DNAPARS` | validated | DNA parsimony, Fitch criterion, stepwise addition + NNI |
| `PROTPARS` | documented difference | Protein parsimony, unordered replacements |
| `DNAPENNY` | partial | DNA parsimony searched exhaustively over addition orders |
| `DNAMOVE` | documented difference | Interactive DNA parsimony |
| `DNACOMP` | validated | DNA compatibility: tree from the largest compatible site set |
| `PARS` | validated | Wagner parsimony for ordered multistate characters |
| `MIX` | validated | Weighted Wagner / Dollo / Camin–Sokal parsimony |
| `PENNY` | partial | Wagner parsimony, intensive search |
| `MOVE` | documented difference | Interactive mixed parsimony |
| `DOLLOP` | validated | Dollo parsimony: one gain, unlimited losses |
| `DOLPENNY` | partial | Dollo parsimony, intensive search |
| `DOLMOVE` | documented difference | Dollo combined with Camin–Sokal |
| `CLIQUE` | validated | Largest clique of mutually compatible characters |

### Likelihood

| Program | Status | What it does |
| --- | --- | --- |
| `DNAML` | validated | DNA ML (JC, F81, K2P, HKY, TN93, GTR; optional Γ rates) |
| `DNAMLK` | documented difference | DNA ML under a strict clock |
| `PROML` | documented difference | Protein ML (Poisson / equal-input) |
| `PROMLK` | documented difference | Protein ML under a strict clock |
| `RESTML` | documented difference | Restriction-site ML (two-state model) |

### Distances and distance trees

| Program | Status | What it does |
| --- | --- | --- |
| `DNADIST` | validated | Pairwise DNA distances (p, JC, K2P, F81, HKY, TN93, GTR, ML) |
| `PROTDIST` | documented difference | Pairwise protein distances |
| `RESTDIST` | validated | Restriction-site distances |
| `GENDIST` | documented difference | Nei, Cavalli-Sforza chord, Reynolds distances |
| `TREEDIST` | documented difference | Robinson–Foulds and branch-score tree distances |
| `NEIGHBOR` | validated | Neighbor joining / UPGMA |
| `FITCH` | validated | Fitch–Margoliash least squares |
| `KITSCH` | validated | Fitch–Margoliash with a clock |
| `CONTML` | documented difference | Brownian-motion restricted ML for continuous data and gene frequencies |
| `CONTRAST` | partial | Felsenstein's independent contrasts |

### Resampling, diagnostics, drawing

| Program | Status | What it does |
| --- | --- | --- |
| `SEQBOOT` | validated | Bootstrap / jackknife replication, optionally searched and summarised |
| `CONSENSE` | validated | Strict, majority-rule or threshold consensus |
| `RETREE` | partial | Reroot, resolve polytomies, rearrange |
| `DNAINVAR` | partial | Cavender's invariants for four taxa |
| `FACTOR` | validated | Multistate to binary recoding through character-state trees |
| `DRAWGRAM` | partial | Rectangular phylogram / cladogram |
| `DRAWTREE` | partial | Radial unrooted phylogram |

## Input formats

| Input | Format |
| --- | --- |
| Sequences | FASTA, or sequential/interleaved PHYLIP (`<taxa> <sites>` first line) |
| Discrete characters | Same, with states `0`–`9` |
| Restriction sites | Same, `1` present and `0` absent |
| Numeric data | Optional `taxa characters` header, then `name value value …` per line |
| Trees | Newick, one per line for multi-tree programs |
| Character-state trees | `#state 3: 0-1 1-2` lines before the sequences (also `0:1 1:2`) |

Missing data (`?`, `-`, `N`, `X`) is skipped when a distance is computed and treated as compatible with any state in parsimony. Phoenix never silently drops an observation: every program reports what it skipped.

## Reading the output

Every result shows a metadata strip: the program and algorithm, its scientific status, a description of the input, the parameters used, the headline statistic, and a timestamp. Exports:

- **Newick** — the tree, with branch lengths.
- **CSV** — a distance matrix.
- **JSON** — the full record: program, status, parameters, input summary, statistics, notes and tree, enough to reproduce the run.
- **Text** — the plain-text report shown on screen.

## Validation

Phoenix is validated against **mathematical references**, not against PHYLIP output. This distinction is deliberate and is stated in [`validation/PHYLIP-DIFFERENTIAL.md`](validation/PHYLIP-DIFFERENTIAL.md).

`npm run validate` runs 62 checks, each recording its expected value, actual value, tolerance and the independent reference used:

- closed-form JC69 and K2P distances;
- exhaustive summation over all 16 internal-state combinations for the likelihood;
- exhaustive enumeration over all internal-state assignments for Dollo parsimony;
- hand-computed Fitch, Wagner and Dollo scores;
- additive distance matrices generated from a known tree;
- simulation on a known topology for the invariant test;
- an independent Gauss–Jordan evaluation of the Brownian-motion restricted likelihood;
- seeded resampling for bootstrap reproducibility.

Results are written to `validation/report.json`. `npm test` adds 88 unit and integration tests, and `npm run browser-check` drives the real interface in Chromium.

## Performance

Measured in a containerised sandbox (see `validation/benchmark.json` for the recorded environment):

| Analysis | Dataset | Time |
| --- | --- | --- |
| `DNADIST` (JC) | 30 taxa × 1000 sites | ~20 ms |
| `NEIGHBOR` (NJ) | 30 taxa | ~24 ms |
| `DNAPARS` | 10 taxa × 500 sites | ~150 ms |
| `DNAPARS` | 30 taxa × 1000 sites | ~7 s |
| `DNAML` (HKY) | 16 taxa × 300 sites | ~15 s |
| `SEQBOOT` + parsimony + `CONSENSE` | 12 taxa × 400 sites, 50 replicates | ~10 s |

Long analyses run in a Web Worker, so the interface stays usable while they compute. Absolute timings will differ on other hardware; treat these as orders of magnitude, not promises.

## Screenshots

Regenerate with `npm run screenshots` (writes to `docs/screenshots/`). Binary images are not committed, so the repository stays text-only and diffable.

## Known limitations

- **No PHYLIP differential comparison has been performed.** Validation is against mathematics. PHYLIP compatibility is not claimed.
- **Searches are heuristic.** `DNAPENNY`, `PENNY` and `DOLPENNY` search exhaustively over addition orders but do not implement branch-and-bound.
- **No empirical amino-acid exchangeability matrices** (JTT, Dayhoff, WAG). Protein likelihood and distances use Poisson or equal-input rates.
- **`CONTML`** maximises the correct restricted likelihood but only relative branch lengths are identifiable; it uses coordinate-wise search rather than Felsenstein's EM algorithm.
- **`DNAINVAR`** implements Cavender's quadratic invariants only, with no significance test, and reports "unresolved" when the data cannot discriminate topologies.
- **No alignment** is performed; sequences must already be aligned.
- Tuned for tens of taxa, not thousands.

## Roadmap

- PHYLIP differential validation where a reference binary is available.
- Empirical exchangeability matrices for protein likelihood.
- Branch-and-bound parsimony.
- Ancestral-state rendering directly on the tree.

## Contributing

Run `npm run check` before opening a pull request. Every bug fix should come with a regression test, and every new program needs a status and a reference note in `phoenix-programs.js` — the compatibility matrix is generated from that registry, so documentation cannot drift.

## License and attribution

Phoenix is an independent reimplementation inspired by PHYLIP, created by Joseph Felsenstein and distributed by the University of Washington. Phoenix bundles no PHYLIP code and claims no compatibility beyond what is demonstrated in `validation/`.

See [`LICENSE`](LICENSE) for the terms of this repository, and [`docs/AUDIT.md`](docs/AUDIT.md), [`docs/METHODS.md`](docs/METHODS.md) and [`docs/COMPATIBILITY.md`](docs/COMPATIBILITY.md) for the scientific record.
