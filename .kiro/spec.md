# Phoenix workspace spec

## Goal
Make the PHYLIP toolchest usable without installing legacy command-line software: a browser-first workspace where sequences become trees, distances and reports.

## Requirements
- Accept pasted or uploaded aligned sequence data and validate it, explaining errors without losing user input.
- Provide the thirty-five PHYLIP programs (parsimony, likelihood, distances, distance trees, resampling, diagnostics, drawing) behind one program browser.
- Each program declares its input type and options; the interface builds the form and the output pane from that declaration.
- Render trees as SVG, offer Newick export, and show a plain-text view of every result.
- Run locally in the browser with no sequence data transmitted to a server.
- Keep the analysis engines free of DOM access and cover them with Node's built-in test runner, so the same files power the browser and the tests.

## Architecture
- `phoenix-core.js` — parsing, character encoding, Newick, tree primitives.
- `phoenix-models.js` — substitution models, transition probabilities, pairwise distances.
- `phoenix-infer.js` — UPGMA, neighbor joining, least squares, parsimony, likelihood.
- `phoenix-support.js` — bootstrap, consensus, tree comparison, ordination, compatibility, drawing.
- `phoenix-programs.js` — the program catalogue and its metadata.
- `programs.js` — the program browser UI; `app.js` remains the quick UPGMA workspace.

## Boundaries
This is not the original PHYLIP software. Searches are heuristic, empirical amino-acid exchangeability matrices are not bundled, CONTML and CONTRAST are Brownian-motion approximations, and alignments remain the user's responsibility. Every program states its own limitations in the interface and in `README.md`. Keep the engines pure and the limitations visible.
