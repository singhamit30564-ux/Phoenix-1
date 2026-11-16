# Pitch and judge preparation

## 30-second pitch

> Phylogenetics is how we work out who is related to whom, and PHYLIP is the
> package that defined the field — a command-line program from 1980. Phoenix is
> a browser reimplementation of thirty-five of its programs: parsimony,
> likelihood, distances, bootstrap, consensus, tree drawing. It runs entirely
> locally, so no sequence data ever leaves your machine. And it carries an
> explicit validation status on every program, because a phylogenetic tree is a
> scientific claim, not a visual effect.

## Two-minute technical explanation

**The problem.** Phylogenetic inference is a numerical discipline with
well-defined algorithms. Most web "phylogeny tools" are either wrappers around
a server or simplified toys. Neither is useful for teaching or for checking
your understanding.

**What we built.** A static, dependency-free application with four engine
modules and one program registry:

- `phoenix-core.js` — parsing, character encoding, Newick, tree primitives.
- `phoenix-models.js` — substitution models and a numerically stable matrix
  exponential (scaling + Taylor + squaring).
- `phoenix-infer.js` — UPGMA, neighbour joining, least squares, parsimony
  (Fitch / Sankoff / Dollo / Camin–Sokal), Felsenstein pruning likelihood.
- `phoenix-support.js` — bootstrap, consensus, invariants, ordination,
  compatibility, drawing.

The registry maps each of the 35 program names onto an implementation and a
scientific status. The UI is generated from that registry, and so is the
compatibility matrix — documentation cannot drift from the code.

**The engineering choices that mattered.**

1. **Likelihood performance.** Recomputing the transition matrix per site per
   edge made protein likelihood take ten seconds. Computing one matrix per edge
   per rate class and vectorising across sites gave the same results roughly ten
   times faster.
2. **Correctness over convenience.** The Dollo gain budget is a *global*
   resource. Optimising each child subtree independently lets two siblings each
   spend a gain, which silently understates the cost. Both implementations had
   that bug; they agreed with each other, which is exactly why the test suite
   now compares against exhaustive enumeration instead.
3. **Responsiveness.** Long analyses run in a Web Worker, so the interface keeps
   responding during a seven-second parsimony run.

**Validation.** Sixty-two checks, each recording an expected value, an actual
value, a tolerance and the reference that produced the expectation — closed
forms, exhaustive enumeration, hand computation, or simulation with a known
generating tree. No check compares an implementation against itself.

## Likely judge questions

**"Is this just PHYLIP recompiled to JavaScript?"**
No. It is an independent reimplementation. No PHYLIP source is bundled, and we
have not compared outputs against PHYLIP — `validation/PHYLIP-DIFFERENTIAL.md`
states that explicitly and documents the protocol we would use.

**"How do you know your trees are right?"**
Two ways. On additivity: neighbour joining reproduces all pairwise path lengths
of an additive matrix exactly (checked to 1e-9). On likelihood: the pruning
recursion is checked against brute-force summation over all 16 internal-state
combinations on four taxa. For parsimony, scores are checked against
hand-computed values and, for Dollo, against exhaustive enumeration of every
internal-state assignment.

**"What about the MOVE/Dollo bug?"**
It was real. The gain budget was optimised per child subtree, so two siblings
could each spend a gain — three gains for the price of one. Fixed by combining
child budgets sequentially; regression tests now compare against exhaustive
enumeration, and a dataset that previously scored 3 (impossible) now scores 4.

**"Which models do you support?"**
Nucleotide: JC69, F81, K2P, HKY, TN93, GTR, with optional discrete-gamma rate
variation. Protein: Poisson and equal-input only — empirical matrices such as
JTT and WAG are not bundled, and the interface says so.

**"How fast is it?"**
Distances on 30 taxa × 1000 sites: ~20 ms. Parsimony on the same data: ~7 s.
DNAML with 16 taxa × 300 sites: ~15 s. All in a browser, on one thread, with the
UI staying responsive.

**"What's missing?"**
Branch-and-bound parsimony, empirical amino-acid matrices, a significance test
for the invariants, and a real PHYLIP differential run. Each is listed in the
compatibility matrix with its program.

**"Why not use a phylogenetics library?"**
The point is the browser as a delivery mechanism for reproducible analysis. A
library would still need a server or a build pipeline; this is one HTML file and
a handful of scripts that work from a USB stick.

**"How do you stop it producing confident nonsense?"**
Every program states its assumptions and its status. Results carry a metadata
strip and export a JSON record of parameters and seed. Bootstrap percentages are
labelled as percentages of replicate trees, not confidence.

## Honest weaknesses to volunteer

- No PHYLIP differential validation.
- Heuristic search only; no optimality guarantee.
- Protein models are limited to two rate schemes.
- Tuned for tens of taxa; thousands would need a different architecture.
