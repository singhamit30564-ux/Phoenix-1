/* phoenix-interactive.js — interactive parsimony sessions (DNAMOVE and MOVE).
   A session holds a working tree, scores it, enumerates the nearest-neighbour
   rearrangements and their score changes, and keeps an undo history. */
(function (root, factory) {
  const api = factory(
    typeof require === 'function' ? require('./phoenix-core.js') : root.PCore,
    typeof require === 'function' ? require('./phoenix-infer.js') : root.PInfer
  );
  if (typeof module === 'object' && module.exports) module.exports = api;
  else root.PInteractive = api;
})(typeof globalThis !== 'undefined' ? globalThis : this, function (C, I) {
  'use strict';

  function findById(node, id) {
    if (node.id === id) return node;
    if (!node.children) return null;
    for (const c of node.children) {
      const r = findById(c, id);
      if (r) return r;
    }
    return null;
  }

  function scorable(tree, names) {
    return C.leaves(tree).length === names.length && names.every(n => C.leaves(tree).some(l => l.name === n));
  }

  function createSession(opts) {
    const { chars, names } = opts;
    const kind = opts.kind || 'dna';
    const weights = opts.weights || { wagner: 1, dollo: 0, caminSokal: 0 };
    const history = [];
    let tree = null;

    function cost(t) {
      if (kind === 'protein') return I.fitchScore(t, chars, names);
      if (weights.dollo || weights.caminSokal) {
        let total = 0;
        const nStates = chars.length ? chars[0][0].length : 10;
        if (weights.wagner) total += weights.wagner * I.sankoffScore(t, chars, names, I.wagnerCost(Math.min(nStates, 10)));
        else total += 0;
        if (weights.dollo) total += weights.dollo * I.dolloScore(t, chars, names);
        if (weights.caminSokal) total += weights.caminSokal * I.sankoffScore(t, chars, names, I.caminSokalCost(Math.min(nStates, 10)));
        return total;
      }
      return I.fitchScore(t, chars, names);
    }

    function reset(candidate) {
      if (candidate) {
        tree = C.cloneTree(candidate);
      } else if (names.length <= 2) {
        tree = { name: null, children: [{ name: names[0], length: 0 }, { name: names[1], length: 0 }] };
      } else {
        tree = I.stepwiseSearch(names, t => -cost(t));
      }
      I.assignIds(tree);
      history.length = 0;
      return tree;
    }

    function currentScore() {
      return cost(tree);
    }

    function moves() {
      const out = [];
      const base = currentScore();
      const internals = [];
      C.postorder(tree, n => { if (n.children && n.children.length === 2) internals.push(n.id); });
      for (const nodeId of internals) {
        for (let which = 0; which < 2; which++) {
          const node = findById(tree, nodeId);
          if (!node) continue;
          const A = node.children[which];
          if (!A || !A.children || A.children.length !== 2) continue;
          for (let sub = 0; sub < 2; sub++) {
            const trial = applyTo(tree, { nodeId, which, sub });
            if (!trial) continue;
            const s = cost(trial);
            const partner = A.children[1 - sub];
            out.push({
              id: `${nodeId}:${which}:${sub}`,
              nodeId, which, sub,
              score: s,
              delta: s - base,
              label: `swap ${A.name || 'clade'} with ${node.children[1 - which].name || 'clade'}`
            });
          }
        }
      }
      return out;
    }

    function applyTo(source, move) {
      const copy = C.cloneTree(source);
      const node = findById(copy, move.nodeId);
      if (!node || !node.children) return null;
      const A = node.children[move.which];
      if (!A || !A.children || A.children.length !== 2) return null;
      const sibling = node.children[1 - move.which];
      const moved = A.children[move.sub];
      A.children[move.sub] = sibling;
      node.children[1 - move.which] = moved;
      return copy;
    }

    function apply(moveId) {
      const move = parseMove(moveId);
      const next = applyTo(tree, move);
      if (!next) return { ok: false, score: currentScore() };
      history.push(C.cloneTree(tree));
      tree = next;
      I.assignIds(tree);
      return { ok: true, score: currentScore() };
    }

    function parseMove(moveId) {
      const [nodeId, which, sub] = String(moveId).split(':').map(Number);
      return { nodeId, which, sub };
    }

    function undo() {
      if (!history.length) return { ok: false, score: currentScore() };
      tree = history.pop();
      I.assignIds(tree);
      return { ok: true, score: currentScore() };
    }

    function greedy() {
      let applied = 0;
      for (let guard = 0; guard < 200; guard++) {
        const improving = moves().filter(m => m.delta < -1e-9).sort((a, b) => a.delta - b.delta)[0];
        if (!improving) break;
        apply(improving.id);
        applied++;
      }
      return { applied, score: currentScore() };
    }

    /** Most-parsimonious reconstruction of each character on the current tree. */
    function siteReconstruction(limit) {
      const idx = new Map(names.map((n, i) => [n, i]));
      const rows = [];
      const max = Math.min(chars.length, limit || 12);
      for (let s = 0; s < max; s++) {
        const col = chars[s];
        const down = new Map();
        let steps = 0;
        C.postorder(tree, node => {
          if (!node.children) {
            const v = col[idx.get(node.name)];
            let mask = 0;
            for (let k = 0; k < v.length; k++) if (v[k] > 0) mask |= (1 << k);
            down.set(node, mask);
            return;
          }
          let m = down.get(node.children[0]);
          for (let k = 1; k < node.children.length; k++) {
            const other = down.get(node.children[k]);
            const both = m & other;
            if (both) m = both;
            else { m |= other; steps++; }
          }
          down.set(node, m);
        });
        rows.push({ site: s + 1, changes: steps, states: describe(down.get(tree)) });
      }
      return rows;
    }

    function describe(mask) {
      const out = [];
      for (let i = 0; i < 16; i++) if (mask & (1 << i)) out.push(i);
      return out;
    }

    reset(opts.initialTree);

    return {
      get tree() { return tree; },
      names,
      kind,
      weights,
      score: currentScore,
      setTree(candidate) { return reset(candidate); },
      reset,
      moves,
      apply,
      undo,
      greedy,
      siteReconstruction,
      canUndo: () => history.length > 0,
      historyLength: () => history.length,
      isScorable: () => scorable(tree, names)
    };
  }

  return { createSession };
});
