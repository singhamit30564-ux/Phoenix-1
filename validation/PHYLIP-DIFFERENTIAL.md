# PHYLIP differential validation

## Status: not performed

No PHYLIP 3.69 executable was available in the environment where Phoenix was
developed, and PHYLIP output is not published in a form that could be used as a
reference for these specific data sets. **No differential comparison against
PHYLIP has been run, and no PHYLIP-equivalence claim is made anywhere in this
repository.**

This file exists so that the gap is explicit and so that the comparison can be
performed later without re-deriving the protocol.

## What Phoenix validation does establish

`validation/validate.js` compares Phoenix against independent mathematical
references: closed-form distance corrections, exhaustive enumeration of internal
states for the likelihood, hand-computed parsimony scores, additive matrices for
neighbor joining, simulations with a known generating topology for the invariant
test, and an independent Gauss–Jordan evaluation of the restricted likelihood.
These checks establish **implementation correctness against mathematical
references**. They do not establish agreement with PHYLIP's implementations,
which differ in optimisers, search strategy, tie-breaking and defaults.

## Protocol for a future differential run

For each program below, run PHYLIP 3.69 and Phoenix on the identical data set
with equivalent options, and record every field.

| Field | Value |
| --- | --- |
| PHYLIP version | 3.69 (record the exact build) |
| Program | e.g. `dnapars` |
| Data set | file name and hash |
| Options | PHYLIP menu settings; Phoenix option values |
| Phoenix output | value under comparison |
| PHYLIP output | value under comparison |
| Tolerance | numeric tolerance used, or "exact" for integer scores/topologies |
| Discrepancy | absolute difference and an explanation |

Programs and the values to compare:

| Program | Compare |
| --- | --- |
| `DNAPARS`, `PARS`, `DOLLOP`, `MIX` | tree length (integer, must match exactly), topology (Robinson–Foulds), tie sets |
| `DNADIST`, `PROTDIST`, `RESTDIST`, `GENDIST` | every matrix entry, with tolerance 1e-4 |
| `DNAML`, `DNAMLK`, `PROML`, `RESTML` | log-likelihood (tolerance 1e-3 for the same topology and similar branch lengths) |
| `NEIGHBOR`, `FITCH`, `KITSCH`, `CONTML` | topology (RF), branch lengths (tolerance 1e-3) |
| `CONSENSE` | retained splits and support percentages (exact) |
| `TREEDIST` | symmetric difference (exact) and branch score (1e-4) |
| `SEQBOOT` | replicate data sets under the same seed are **not** expected to match; compare the distribution of support values instead |
| `DNAINVAR` | invariant values per topology |
| `FACTOR` | recoded matrix and the character mapping |

Expected sources of legitimate difference, to be recorded rather than treated as
bugs:

- Heuristic search order: PHYLIP and Phoenix may return different equally
  parsimonious trees.
- Log-likelihood optimisers: different convergence and branch-length tolerance.
- `CONTML`: PHYLIP uses an EM algorithm, Phoenix coordinate-wise search, on the
  same restricted likelihood.
- Defaults: for example PHYLIP's transition/transversion ratio and its
  amino-acid step matrix.
- Randomisation: bootstrap and jumble seeds do not correspond between programs.

## Rule

Do not add a row to this table without the recorded PHYLIP version, the data
set, the inputs and the outputs. A claim of compatibility without that evidence
would be worse than the current explicit statement that it has not been tested.
