# Scientific methods

This document states what each class of calculation in Phoenix actually
computes, the assumptions behind it, and its limitations. Each program also
carries a one-line status and reference in the interface
(`docs/COMPATIBILITY.md` is generated from those).

## Distance correction

For two aligned sequences, let *p* be the fraction of comparable sites that
differ. Sites that are ambiguous or missing in either sequence are skipped for
that pair.

- **p-distance** — *d* = *p*.
- **Jukes–Cantor** (Jukes & Cantor 1969) — *d* = −(3/4) ln(1 − 4*p*/3), defined
  only for *p* < 0.75.
- **Kimura 2-parameter** (Kimura 1980) — with *P* transition and *Q*
  transversion proportions, *d* = ½ ln(1/(1−2*P*−*Q*)) + ¼ ln(1/(1−2*Q*)).
- **F81, HKY, TN93, GTR, and maximum likelihood** — the distance is the value of
  *t* that maximises the two-taxon likelihood under the corresponding rate
  matrix. For F81/HKY/TN93/GTR this is the maximum-likelihood distance, not a
  closed-form correction; when the model coincides with JC or K2P the numerical
  answer agrees with the closed form, which is how it is validated.

## Substitution models

Rate matrices are built as *Q*<sub>ij</sub> = *r*<sub>ij</sub>π<sub>j</sub> with
the diagonal set so that rows sum to zero, then scaled so that the expected rate
−Σπ<sub>i</sub>*Q*<sub>ii</sub> equals one. Branch lengths are therefore
expected substitutions per site.

Transition probabilities come from *P*(*t*) = exp(*Qt*), computed by scaling and
squaring with a Taylor series. The implementation is checked against the
semigroup property *P*(*t*₁)*P*(*t*₂) = *P*(*t*₁+*t*₂) and against row sums of
one.

Rate variation across sites uses a discrete-gamma approximation with mean 1.

## Likelihood

Felsenstein's pruning algorithm gives the likelihood of an alignment on a tree.
For each site, partial likelihoods are computed postorder; ambiguous or missing
states contribute a vector of ones (they are unobserved, not absent).

Branch lengths are optimised coordinate-wise by golden-section maximisation of
the log-likelihood. Topologies are refined by nearest-neighbour interchange,
accepting a rearrangement only if the likelihood increases.

The clock variants reparameterise the tree by node heights and optimise the
heights under the constraint that no parent is younger than its children, which
produces an ultrametric tree.

Constrained-optimisation results depend on the starting tree. A different
addition order can reach a different local optimum; this is reported as a
heuristic search, never as a guarantee.

## Parsimony

- **Fitch** (unordered): postorder set intersection; a union costs one step.
- **Sankoff** with an arbitrary step matrix: handles Wagner (ordered, cost
  |*i*−*j*|) and Camin–Sokal (irreversible, 0 → 1 only).
- **Dollo**: dynamic programming over (gains used ≤ 1, node state), with each
  loss costing one step. States above zero are read as "present".
- **Mixed (MIX)**: a user-weighted sum of the three criteria.

Searches use stepwise addition followed by nearest-neighbour interchange and
random restarts. `DNAPENNY`, `PENNY` and `DOLPENNY` search exhaustively over
addition orders; they do not implement branch-and-bound with bounds.

## Character compatibility

Two characters are compatible when the state sets induced by one character on
the other form a laminar family (Estabrook & McMorris). `CLIQUE` finds a maximum
clique of mutually compatible characters by Bron–Kerbosch; `DNACOMP` does the
same for sites and builds the tree implied by the compatible set.

## Tree building from distances

- **UPGMA** joins the closest clusters using size-weighted average distances and
  assumes a clock, producing an ultrametric tree.
- **Neighbor joining** (Saitou & Nei 1987) is the minimum-evolution criterion
  with the *Q* matrix correction; it recovers an additive matrix exactly, which
  is how it is validated.
- **Fitch–Margoliash** minimises Σ*w*<sub>ij</sub>(*d*<sub>ij</sub> −
  *t*<sub>ij</sub>)², with *w* = 1/*d*² by default. Branch lengths are refit in
  closed form and the topology is refined by NNI.
- **KITSCH** is the same criterion under a clock, optimising node heights.

## Brownian-motion restricted likelihood (CONTML)

Continuous characters and gene frequencies are modelled as Brownian motion on
the tree (Cavalli-Sforza & Edwards 1967; Felsenstein 1973, 1981c). Tip vectors
are multivariate normal with covariance σ²*C*, where *C*<sub>ij</sub> is the
distance from the root to the most recent common ancestor of *i* and *j*.

Both the ancestral state and σ² are profiled out, giving the restricted (REML)
likelihood

  −2 ln L = *nd* ln 2π + *d* ln|*C*| + *d* ln(1ᵀ*C*⁻¹1) + dof·ln(RSS/dof) + dof

The ln(1ᵀ*C*⁻¹1) term is what makes the value invariant to where the tree is
rooted; the validation suite checks both root- and scale-invariance numerically
and cross-checks the value against an independent Gauss–Jordan evaluation.

Because the likelihood is invariant to rescaling all branch lengths while
rescaling σ² inversely, only **relative** branch lengths are identifiable.
Phoenix reports lengths in units of the rate. Gene frequencies are transformed
with arcsin √*p* first, the variance-stabilising transformation of
Cavalli-Sforza & Edwards.

## Invariants (DNAINVAR)

For four taxa, reduce nucleotides to purines and pyrimidines and form the 4×4
matrix *M* whose rows are the state pairs at taxa 1 and 2 and whose columns are
the state pairs at taxa 3 and 4. For the topology 12|34 this matrix factorises
as *U P V* with *U* of size 4×2 and *V* of size 2×4, so it has rank at most two
and every 3×3 minor vanishes in expectation for any branch lengths (Cavender &
Felsenstein 1987). The other two topologies impose the same condition on their
own pairings. Phoenix reports the invariant residual — the Euclidean norm of the
sixteen 3×3 minors — for each topology. This is Cavender's quadratic invariant,
not Lake's linear invariant, and no significance test is applied.

## Recoding (FACTOR)

Each multistate character is given a character-state tree (a linear chain of the
observed states unless one is supplied with `#state N: a-b b-c` lines). Every
edge of that tree becomes one binary character, set to 1 for taxa on the derived
side. Since a state tree with *k* states has *k*−1 edges and every state
corresponds to a unique pattern of edges, the recoding is exactly reversible, and
the implementation verifies the round trip for every observed state.

## Support

Bootstrap and jackknife resample alignment columns (not taxa) with a seeded
random-number generator, so a fixed seed reproduces the same replicates.
Consensus trees are built by counting splits, keeping those above a threshold,
and adding them greedily while they remain mutually compatible (the
majority-rule or extended consensus). Internal labels are percentages of input
trees containing that split — they describe the data, not a posterior
probability of the topology.

## References

- Cavender, J. A. & Felsenstein, J. (1987) *Mol. Biol. Evol.* 4:625–640.
- Cavalli-Sforza, L. L. & Edwards, A. W. F. (1967) *Evolution* 21:550–570.
- Felsenstein, J. (1973) *Genetics* 79:155–170.
- Felsenstein, J. (1981) *J. Mol. Evol.* 17:368–376.
- Felsenstein, J. (1981c) *Evolution* 35:873–890 (restricted likelihood / EM).
- Fitch, W. M. (1971) *Syst. Zool.* 20:406–416.
- Jukes, T. H. & Cantor, C. R. (1969) in *Mammalian Protein Metabolism*.
- Kimura, M. (1980) *J. Mol. Evol.* 16:111–120.
- Saitou, N. & Nei, M. (1987) *Mol. Biol. Evol.* 4:406–425.
- Sankoff, D. (1975) *Syst. Zool.* 24:424–435.
- Robinson, D. F. & Foulds, L. R. (1981) *Math. Biosci.* 53:131–147.
