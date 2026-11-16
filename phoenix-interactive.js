/* phoenix-interactive.js — interactive parsimony sessions (DNAMOVE and MOVE).
 *
 * A session holds a working tree, scores it under an explicitly selected
 * criterion, enumerates nearest-neighbour rearrangements with their score
 * changes, and keeps an undo history.
 *
 * Scoring is dispatched on the criterion, never on a fall-through:
 *
 *   DNA / protein  -> Fitch (unordered), as DNAMOVE and PROTPARS do
 *   discrete data  -> the weighted Wagner / Dollo / Camin-Sokal objective,
 *                     which is what MOVE (and MIX) use
 *
 * For discrete data the score and the ancestral reconstruction come from one
 * dynamic program over (gains used, node state), so the states shown are the
 * states that actually achieve the reported score.
 */
(function (root, factory) {
  const api = factory(
    typeof require === 'function' ? require('./phoenix-core.js') : root.PCore,
    typeof require === 'function' ? require('./phoenix-infer.js') : root.PInfer
  );
  if (typeof module === 'object' && module.exports) module.exports = api;
  else root.PInteractive = api;
})(typeof globalThis !== 'undefined' ? globalThis : this, function (C, I) {
  'use strict';

  const INF = 1e12;
  const G = 2;            // gain budget dimension: 0 or 1 gain used
  const MAX_STATES = 10;  // discrete states are 0..9

  function findById(node, id) {
    if (node.id === id) return node;
    if (!node.children) return null;
    for (const c of node.children) {
      const r = findById(c, id);
      if (r) return r;
    }
    return null;
  }

  /**
   * Resolve the requested weights into a criterion. Throws rather than
   * silently substituting a different criterion when the configuration is
   * ambiguous or empty.
   */
  function resolveCriterion(kind, rawWeights) {
    if (kind === 'dna' || kind === 'protein') {
      return { type: 'fitch' };
    }
    // Only a completely absent weight specification falls back to Wagner.
    // An explicit { dollo: 1 } means Dollo only, never "Dollo plus Wagner".
    const w = rawWeights || { wagner: 1, dollo: 0, caminSokal: 0 };
    const weights = {
      wagner: Number(w.wagner || 0),
      dollo: Number(w.dollo || 0),
      caminSokal: Number(w.caminSokal || 0)
    };
    for (const [name, value] of Object.entries(weights)) {
      if (!Number.isFinite(value) || value < 0) {
        throw new Error(`The ${name} weight must be a number greater than or equal to zero.`);
      }
    }
    const total = weights.wagner + weights.dollo + weights.caminSokal;
    if (!(total > 0)) {
      throw new Error('At least one of the Wagner, Dollo or Camin–Sokal weights must be greater than zero.');
    }
    const active = Object.entries(weights).filter(([, v]) => v > 0).map(([k]) => k);
    if (active.length === 1) return { type: active[0], weights };
    return { type: 'mixed', weights };
  }

  /**
   * Cost of one state transition on one branch. Each active criterion
   * contributes its own cost; Camin-Sokal makes downward moves impossible.
   */
  function edgeCost(from, to, weights) {
    let cost = 0;
    if (weights.wagner) cost += weights.wagner * Math.abs(from - to);
    if (weights.caminSokal) {
      if (to < from) return INF;          // irreversible: no state may revert
      cost += weights.caminSokal * (to - from);
    }
    if (weights.dollo) {
      const gain = from === 0 && to > 0 ? 1 : 0;
      const loss = from > 0 && to === 0 ? 1 : 0;
      cost += weights.dollo * (gain + loss);
    }
    return cost;
  }

  function isGain(from, to) {
    return from === 0 && to > 0 ? 1 : 0;
  }

  /**
   * One dynamic program over (gains used, state) that serves every criterion:
   * with only Wagner active it reduces to Sankoff with ordered costs, with only
   * Dollo active to the binary presence/absence Dollo criterion, and with a
   * mixture to the weighted sum. Returns the cost and a MAP reconstruction.
   */
  function solveCharacter(tree, column, names, weights, nStates) {
    const index = new Map(names.map((n, i) => [n, i]));
    // The gain budget is a global resource and only exists when Dollo is
    // active. Tracking it for Wagner-only or Camin-Sokal-only runs would
    // wrongly forbid solutions that use several gains.
    const trackGains = weights.dollo > 0;
    const dims = trackGains ? 2 : 1;
    const size = dims * nStates;
    const dp = new Map();
    const plan = new Map();

    C.postorder(tree, node => {
      const table = new Float64Array(size).fill(INF);
      if (!node.children) {
        const v = column[index.get(node.name)];
        const allowed = [];
        for (let s = 0; s < nStates; s++) if (v[s] > 0) allowed.push(s);
        if (!allowed.length) table.fill(0);          // unobserved: any state, no gain
        else for (const s of allowed) table[s] = 0;  // a tip never spends a gain itself
        dp.set(node, table);
        plan.set(node, null);
        return;
      }

      const plans = new Array(size).fill(null);
      for (let g = 0; g < dims; g++) {
        for (let s = 0; s < nStates; s++) {
          // Children are combined one after another so that a gain spent in
          // one subtree is no longer available in its siblings.
          let acc = new Float64Array(dims).fill(INF);
          acc[0] = 0;
          const perChild = [];
          let feasible = true;
          for (const child of node.children) {
            const ct = dp.get(child);
            const next = new Float64Array(dims).fill(INF);
            const decision = new Array(dims).fill(null);
            for (let used = 0; used <= g; used++) {
              if (acc[used] >= INF) continue;
              for (let c = 0; c < nStates; c++) {
                const step = edgeCost(s, c, weights);
                if (step >= INF) continue;
                // Only Dollo spends budget on a 0 -> derived transition; for the
                // other criteria an upward move is just an ordinary step.
                const gain = trackGains ? isGain(s, c) : 0;
                for (let gc = 0; used + gain + gc <= g; gc++) {
                  const v = ct[gc * nStates + c];
                  if (v >= INF) continue;
                  const candidate = acc[used] + step + v;
                  const nu = used + gain + gc;
                  if (candidate < next[nu] - 1e-12) {
                    next[nu] = candidate;
                    decision[nu] = { c, budget: gc, from: used };
                  }
                }
              }
            }
            if (next.every(x => x >= INF)) { feasible = false; break; }
            perChild.push(decision);
            acc = next;
          }
          if (!feasible) continue;
          let bestUsed = -1, bestCost = INF;
          for (let used = 0; used <= g; used++) {
            if (acc[used] < bestCost - 1e-12) { bestCost = acc[used]; bestUsed = used; }
          }
          if (bestUsed < 0) continue;
          table[g * nStates + s] = bestCost;
          plans[g * nStates + s] = { used: bestUsed, perChild };
        }
      }
      dp.set(node, table);
      plan.set(node, plans);
    });

    const rootTable = dp.get(tree);
    let bestCost = INF, bestG = -1, bestS = -1;
    for (let g = 0; g < dims; g++) {
      for (let s = 0; s < nStates; s++) {
        const v = rootTable[g * nStates + s];
        if (v < bestCost - 1e-12) { bestCost = v; bestG = g; bestS = s; }
      }
    }
    if (bestG < 0 || !isFinite(bestCost) || bestCost >= INF) {
      return { cost: INF, states: new Map(), rootState: -1, changes: 0, gains: 0, losses: 0 };
    }

    const states = new Map();
    const edges = [];
    (function trace(node, g, s) {
      states.set(node, s);
      if (!node.children) return;
      const entry = plan.get(node)[g * nStates + s];
      if (!entry) return;
      let used = entry.used;
      // Walk the per-child decisions backwards: each records the usage it
      // inherited and the usage it produced.
      for (let k = node.children.length - 1; k >= 0; k--) {
        const decision = entry.perChild[k][used];
        if (!decision) break;
        const child = node.children[k];
        edges.push([s, decision.c]);
        trace(child, decision.budget, decision.c);
        used = decision.from;
      }
    })(tree, bestG, bestS);

    let gains = 0, losses = 0;
    for (const [from, to] of edges) {
      if (from === 0 && to > 0) gains++;
      else if (from > 0 && to === 0) losses++;
    }

    return {
      cost: bestCost,
      states,
      rootState: bestS,
      changes: edges.filter(([from, to]) => from !== to).length,
      gains,
      losses
    };
  }

  // ------------------------------------------------------------- Fitch path

  function fitchSolve(tree, column, names) {
    const index = new Map(names.map((n, i) => [n, i]));
    const down = new Map();
    const finalSet = new Map();
    let steps = 0;
    C.postorder(tree, node => {
      if (!node.children) {
        const v = column[index.get(node.name)];
        let mask = 0;
        for (let s = 0; s < v.length; s++) if (v[s] > 0) mask |= (1 << s);
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
    let changes = 0;
    (function pass(node, parentMask) {
      const own = down.get(node);
      const chosen = parentMask !== undefined && (own & parentMask) ? (own & parentMask) : own;
      finalSet.set(node, chosen);
      if (node.children) {
        node.children.forEach(child => {
          const childMask = down.get(child);
          if (!(childMask & chosen)) changes++;
          pass(child, chosen);
        });
      }
    })(tree, undefined);
    return { steps, states: finalSet, changes };
  }

  function maskToStates(mask) {
    const out = [];
    for (let i = 0; i < 16; i++) if (mask & (1 << i)) out.push(i);
    return out;
  }

  // ----------------------------------------------------------- the session

  function createSession(opts) {
    const { chars, names } = opts;
    const kind = opts.kind || 'dna';
    const criterion = resolveCriterion(kind, opts.weights);
    const nStates = MAX_STATES;
    const history = [];
    let tree = null;

    function costOf(t) {
      if (criterion.type === 'fitch') return I.fitchScore(t, chars, names);
      let total = 0;
      for (const col of chars) {
        total += solveCharacter(t, col, names, criterion.weights, nStates).cost;
      }
      return total;
    }

    function reset(candidate) {
      if (candidate) {
        tree = C.cloneTree(candidate);
      } else if (names.length <= 2) {
        tree = { name: null, children: [{ name: names[0], length: 0 }, { name: names[1], length: 0 }] };
      } else {
        tree = I.stepwiseSearch(names, t => -costOf(t));
      }
      I.assignIds(tree);
      history.length = 0;
      return tree;
    }

    function currentScore() {
      return costOf(tree);
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
            const score = costOf(trial);
            out.push({
              id: `${nodeId}:${which}:${sub}`,
              nodeId, which, sub,
              score,
              delta: score - base,
              label: `swap ${A.name || 'clade'} with ${node.children[1 - which].name || 'clade'}`
            });
          }
        }
      }
      return out;
    }

    function apply(moveId) {
      const [nodeId, which, sub] = String(moveId).split(':').map(Number);
      const next = applyTo(tree, { nodeId, which, sub });
      if (!next) return { ok: false, score: currentScore() };
      history.push(C.cloneTree(tree));
      tree = next;
      I.assignIds(tree);
      return { ok: true, score: currentScore() };
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

    /**
     * One most-parsimonious reconstruction per character on the current tree.
     * Runs the same criterion as the score, so `cost` here sums to the score.
     * Ties between equally parsimonious reconstructions are broken toward the
     * lowest state index.
     */
    function siteReconstruction(limit) {
      const max = Math.min(chars.length, limit || 12);
      const rows = [];
      for (let s = 0; s < max; s++) {
        const col = chars[s];
        if (criterion.type === 'fitch') {
          const solved = fitchSolve(tree, col, names);
          rows.push({
            site: s + 1,
            cost: solved.steps,
            changes: solved.changes,
            rootStates: maskToStates(solved.states.get(tree))
          });
        } else {
          const solved = solveCharacter(tree, col, names, criterion.weights, nStates);
          rows.push({
            site: s + 1,
            cost: Number(solved.cost.toFixed(3)),
            changes: solved.changes,
            gains: solved.gains,
            losses: solved.losses,
            rootStates: solved.rootState >= 0 ? [solved.rootState] : []
          });
        }
      }
      return rows;
    }

    function criterionLabel() {
      if (criterion.type === 'fitch') return 'Fitch (unordered)';
      if (criterion.type === 'mixed') {
        const w = criterion.weights;
        return `weighted Wagner ${w.wagner} / Dollo ${w.dollo} / Camin–Sokal ${w.caminSokal}`;
      }
      return criterion.type === 'caminSokal'
        ? 'Camin–Sokal (irreversible)'
        : criterion.type.charAt(0).toUpperCase() + criterion.type.slice(1);
    }

    reset(opts.initialTree);

    return {
      get tree() { return tree; },
      names,
      kind,
      criterion: criterion.type,
      criterionLabel: criterionLabel(),
      score: currentScore,
      setTree(candidate) { return reset(candidate); },
      reset,
      moves,
      apply,
      undo,
      greedy,
      siteReconstruction,
      canUndo: () => history.length > 0,
      historyLength: () => history.length
    };
  }

  return { createSession, resolveCriterion, edgeCost, solveCharacter, fitchSolve };
});
