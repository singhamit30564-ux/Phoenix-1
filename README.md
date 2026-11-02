# Phoenix

A browser-based, PHYLIP-inspired phylogenetics workspace. Paste or upload aligned DNA sequences in FASTA format, select a distance model, build a UPGMA tree, and export it as Newick. Everything runs locally in your browser; no sequences are sent to a server.

## Run

Open `index.html` in a modern browser, or serve the folder with `python3 -m http.server 8000` and visit `http://localhost:8000`. No build step or dependencies are required. The interface uses Google Fonts when available, with local font fallbacks.

## Analysis and limitations

- Input must contain 2–100 uniquely named, equal-length, **already aligned** DNA sequences. Supported characters are A, C, G, T, N, `?`, and `-`; missing or ambiguous sites are excluded pairwise. Pairs with no comparable sites are rejected.
- Choose raw p-distance (fraction of comparable sites that differ) or Jukes–Cantor correction. The latter is undefined at or above 75% divergence and is rejected in that case.
- UPGMA produces a rooted, distance-based tree using cluster-size-weighted average distances. It assumes a molecular clock; it is not a substitute for maximum-likelihood inference or a full PHYLIP installation. No alignment, bootstrap support, or statistical confidence estimates are provided. Results are intended for exploration, not clinical use.
- Exported Newick branch lengths represent distance units. Input names are quoted and escaped. The visualization abbreviates long names, while export preserves them.

## Tests

Run `node --test tests/*.test.js` with Node.js 18 or newer.

## Credits

Inspired by [PHYLIP](https://phylipweb.github.io/phylip/), created by Joseph Felsenstein. Phoenix is an independent reimplementation of a small subset of phylogenetic functionality; it does not bundle PHYLIP or claim compatibility with all its programs.
