# Phoenix workspace spec

## Goal
Make a small, approachable phylogenetics workflow usable without installing legacy command-line software.

## Requirements
- Accept pasted or uploaded aligned DNA FASTA and provide an example dataset.
- Validate input and explain errors without losing user data.
- Compute pairwise p-distance or Jukes–Cantor distances, then infer a UPGMA tree.
- Display a rooted tree and export the result in Newick format.
- Run locally in the browser with no sequence data transmitted to a server.

## Boundaries
This is not the original PHYLIP software or a comprehensive replacement. UPGMA assumes a molecular clock; alignments, statistical support, and more advanced inference are outside the initial scope. Keep the analysis engine separate from the interface and test it with Node's built-in test runner.
