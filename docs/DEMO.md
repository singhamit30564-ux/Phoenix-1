# Five-minute live demonstration

A script for presenting Phoenix. Timings assume a laptop and a projector.

## Before you start

```bash
python3 -m http.server 8000
# open http://localhost:8000 in a browser at ~1280px wide
```

Have a second tab open on `docs/COMPATIBILITY.md` if a judge asks about a
specific program's status.

## 0:00 — The problem (30 seconds)

> "PHYLIP is the package that taught a generation of biologists phylogenetics.
> It's also a command-line program from 1980 with a menu you have to memorise.
> Phoenix puts the same methods in a browser tab — and, unlike most web
> reimplementations, it tells you exactly how far each one has been validated."

Point at the header: **35 programs**, all running locally.

## 0:30 — Demo A: sequences to tree (60 seconds)

1. Scroll to **GUIDED DEMOS**, click **Run demo** on card A.
2. While it runs: "Five primate sequences, 32 sites. Two steps: pairwise
   distances, then neighbour joining."
3. Point at the metadata strip: program, status, input summary, parameters,
   timestamp.
4. Point at the distance matrix in the text pane: "Jukes–Cantor correction —
   the diagonal is zero, and the pairs are symmetric, which is a cheap sanity
   check on the implementation."
5. Click **↓ Newick** to show the export.

## 1:30 — Demo B: parsimony and ancestral states (60 seconds)

1. Click **Run demo** on card B.
2. "Wagner parsimony: states are ordered, so a change of two states costs two
   steps. The score is the total number of steps."
3. Scroll to the character table: "For each character we get the cost, how many
   branches change, and the reconstructed state at the root. The costs sum to
   the score, because the score and the reconstruction come from the same
   dynamic program."
4. Optional: go to the **MOVE** program in the sidebar, run it, and apply a
   rearrangement. "Interactive parsimony — click a rearrangement, watch the
   score change, undo it."

## 2:30 — Demo C: bootstrap support (60 seconds)

1. Click **Run demo** on card C.
2. "Fifty bootstrap replicates, each searched with parsimony, summarised as a
   majority-rule consensus."
3. "The labels are percentages of replicate trees containing that split. All
   three clades sit at 100% — and they should, because this alignment was
   *simulated* on a known tree. That's the point: we can check the answer."
4. Open the JSON export to show the seed and parameters recorded for
   reproducibility.

## 3:30 — Why it's trustworthy (60 seconds)

> "Every program carries a status. Fifteen are validated — meaning a number they
> produce has been checked against something independent: a closed-form formula,
> an exhaustive enumeration, or a simulation with a known answer. Twelve have
> documented differences from PHYLIP. Eight are partial, and the matrix says
> exactly what's missing."

If asked for evidence: run `npm run validate` (62 checks) or open
`validation/report.json`.

## 4:30 — Close (30 seconds)

> "It's a static site — no build step, no dependencies, no server. Clone it and
> serve the folder. The scientific claims are in the repository, not just in the
> slides."

---

## Failure drill

If something does not render during the demo:

- Reload the page; demos are deterministic and the seed is fixed.
- If the page was opened as `file://`, the Web Worker is blocked (browsers
  forbid it); analyses still run on the main thread, just without the
  responsiveness benefit. Serve over HTTP to avoid this.
- Demo C takes up to ~10 seconds; the status line shows elapsed time. Do not
  click twice — the run is not cancellable.

## What not to say

- Do not claim Phoenix reproduces PHYLIP output. It has not been compared
  against PHYLIP.
- Do not claim the heuristics find optimal trees; they find good ones.
- Do not claim the bootstrap percentages are confidence values.
