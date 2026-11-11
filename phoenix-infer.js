/* phoenix-infer.js — tree inference: UPGMA, neighbor joining, least squares
   (FITCH/KITSCH), parsimony (Wagner / Dollo / Camin-Sokal) and maximum likelihood. */
(function (root, factory) {
  const api = factory(
    typeof require === 'function' ? require('./phoenix-core.js') : root.PCore,
    typeof require === 'function' ? require('./phoenix-models.js') : root.PModels
  );
  if (typeof module === 'object' && module.exports) module.exports = api;
  else root.PInfer = api;
})(typeof globalThis !== 'undefined' ? globalThis : this, function (PCore, PModels) {
  'use strict';

  const { leaves, cloneTree, parentMap, postorder } = PCore;
  const INF = 1e9;

  // ------------------------------------------------------------- clustering

  function upgma(D, names) {
    let clusters = names.map((n, i) => ({ name: n, size: 1, height: 0, idx: i }));
    let d = D.map(r => r.slice());
    while (clusters.length > 1) {
      let x = 0, y = 1, best = Infinity;
      for (let i = 0; i < clusters.length; i++) {
        for (let j = i + 1; j < clusters.length; j++) {
          if (d[i][j] < best) { best = d[i][j]; x = i; y = j; }
        }
      }
      if (!isFinite(best)) throw new Error('Distance matrix has unreachable values.');
      const a = clusters[x], b = clusters[y];
      const merged = { children: [a, b], height: best / 2, size: a.size + b.size, name: null };
      const keep = clusters.map((_, i) => i).filter(i => i !== x && i !== y);
      const nextD = keep.map(i => keep.map(j => d[i][j]));
      const averages = keep.map(i => (d[x][i] * a.size + d[y][i] * b.size) / merged.size);
      nextD.forEach((row, i) => row.push(averages[i]));
      nextD.push([...averages, 0]);
      clusters = [...keep.map(i => clusters[i]), merged];
      d = nextD;
    }
    const root = clusters[0];
    (function setLengths(n) {
      if (!n.children) return;
      n.length = n.length === undefined ? null : n.length;
      for (const c of n.children) { c.length = Math.max(0, n.height - (c.height || 0)); setLengths(c); }
    })(root);
    root.length = null;
    return root;
  }

  function neighborJoining(D, names) {
    const n = names.length;
    if (n === 2) {
      return { name: null, children: [{ name: names[0], length: D[0][1] / 2 }, { name: names[1], length: D[0][1] / 2 }] };
    }
    let nodes = names.map(nm => ({ name: nm, children: null, length: null }));
    let d = D.map(r => r.slice());
    while (nodes.length > 2) {
      const m = nodes.length;
      const r = new Array(m).fill(0);
      for (let i = 0; i < m; i++) for (let j = 0; j < m; j++) if (i !== j) r[i] += d[i][j];
      let bx = 0, by = 1, best = Infinity;
      for (let i = 0; i < m; i++) {
        for (let j = i + 1; j < m; j++) {
          const q = (m - 2) * d[i][j] - r[i] - r[j];
          if (q < best) { best = q; bx = i; by = j; }
        }
      }
      const li = 0.5 * d[bx][by] + (r[bx] - r[by]) / (2 * (m - 2));
      const lj = d[bx][by] - li;
      nodes[bx].length = Math.max(0, li);
      nodes[by].length = Math.max(0, lj);
      const parent = { name: null, children: [nodes[bx], nodes[by]], length: null };
      const newD = [];
      for (let k = 0; k < m; k++) {
        if (k === bx || k === by) continue;
        newD.push((d[bx][k] + d[by][k] - d[bx][by]) / 2);
      }
      const reduced = [];
      const keepIdx = [];
      for (let k = 0; k < m; k++) if (k !== bx && k !== by) keepIdx.push(k);
      const mat = keepIdx.map((a, i) => keepIdx.map((b, j) => (i === j ? 0 : d[a][b])));
      mat.forEach((row, i) => row.push(newD[i]));
      mat.push([...newD, 0]);
      reduced.push(mat);
      const nextNodes = keepIdx.map(k => nodes[k]);
      nextNodes.push(parent);
      nodes = nextNodes;
      d = reduced[0];
    }
    const root = { name: null, children: [nodes[0], nodes[1]], length: null };
    const half = (d[0][1] || 0) / 2;
    nodes[0].length = Math.max(0, half);
    nodes[1].length = Math.max(0, half);
    return root;
  }

  function assignIds(node) {
    let i = 0;
    (function walk(n) { n.id = i++; if (n.children) n.children.forEach(walk); })(node);
    return node;
  }

  function findById(node, id) {
    if (node.id === id) return node;
    if (!node.children) return null;
    for (const c of node.children) { const r = findById(c, id); if (r) return r; }
    return null;
  }

  // ----------------------------------------------------------- path lengths

  function treeIndices(root) {
    const names = leaves(root).map(l => l.name);
    const leafByName = new Map(leaves(root).map(l => [l.name, l]));
    const parent = parentMap(root);
    const pathToRoot = new Map();
    const pathOf = n => {
      if (pathToRoot.has(n)) return pathToRoot.get(n);
      const out = n === root ? [root] : [...pathOf(parent.get(n)), n];
      pathToRoot.set(n, out);
      return out;
    };
    return { names, leafByName, parent, pathOf };
  }

  function mrca(a, b, ctx) {
    const pa = ctx.pathOf(a), pb = ctx.pathOf(b);
    let i = 0;
    while (i < pa.length && i < pb.length && pa[i] === pb[i]) i++;
    return pa[i - 1];
  }

  function pathLength(a, b, ctx) {
    const pa = ctx.pathOf(a), pb = ctx.pathOf(b);
    const setB = new Set(pb);
    const shared = new Set(pa.filter(n => setB.has(n)));
    let total = 0;
    for (const n of pa) if (!shared.has(n)) total += n.length || 0;
    for (const n of pb) if (!shared.has(n)) total += n.length || 0;
    return total;
  }

  // ---------------------------------------------- least-squares tree fitting

  function leastSquaresFit(D, names, opts = {}) {
    const clock = !!opts.clock;
    const wt = opts.weights === 'equal' ? () => 1
      : opts.weights === '1/d' ? d => (d > 0 ? 1 / d : 1)
        : d => (d > 0 ? 1 / (d * d) : 1);
    let tree = clock ? upgma(D, names) : neighborJoining(D, names);
    assignIds(tree);
    const ctx0 = treeIndices(tree);
    void ctx0;
    for (let sweep = 0; sweep < 15; sweep++) {
      if (clock) optimizeHeights(tree, D, wt);
      else optimizeBranchLengthsLS(tree, D, wt);
      const improved = nniImprove(tree, t => {
        const c = treeIndices(t);
        return -weightedSS(t, D, wt, c);
      });
      if (!improved) break;
    }
    if (clock) optimizeHeights(tree, D, wt);
    else optimizeBranchLengthsLS(tree, D, wt);
    const ctx = treeIndices(tree);
    return { tree, ss: weightedSS(tree, D, wt, ctx) };
  }

  function weightedSS(tree, D, wt, ctx) {
    const { names, leafByName } = ctx;
    let ss = 0;
    for (let i = 0; i < names.length; i++) {
      for (let j = i + 1; j < names.length; j++) {
        const w = wt(D[i][j]);
        const diff = D[i][j] - pathLength(leafByName.get(names[i]), leafByName.get(names[j]), ctx);
        ss += w * diff * diff;
      }
    }
    return ss;
  }

  function optimizeBranchLengthsLS(tree, D, wt) {
    for (let round = 0; round < 5; round++) {
      const ctx = treeIndices(tree);
      const { names, leafByName, parent } = ctx;
      for (const [child] of parent) {
        const below = new Set(leaves(child).map(l => l.name));
        let num = 0, den = 0;
        child.length = 0;
        for (let i = 0; i < names.length; i++) {
          for (let j = i + 1; j < names.length; j++) {
            if (below.has(names[i]) === below.has(names[j])) continue;
            const w = wt(D[i][j]);
            const implied = pathLength(leafByName.get(names[i]), leafByName.get(names[j]), ctx);
            num += w * (D[i][j] - implied);
            den += w;
          }
        }
        child.length = den > 0 ? Math.max(0, num / den) : 0;
      }
    }
  }

  function optimizeHeights(tree, D, wt) {
    const ctx = treeIndices(tree);
    const { names, leafByName } = ctx;
    const heights = new Map();
    postorder(tree, n => heights.set(n, n.children ? Math.max(...n.children.map(c => heights.get(c) + (c.length || 0))) : 0));
    const internals = [];
    postorder(tree, n => { if (n.children) internals.push(n); });
    const pairs = [];
    for (let i = 0; i < names.length; i++) {
      for (let j = i + 1; j < names.length; j++) {
        pairs.push([i, j, mrca(leafByName.get(names[i]), leafByName.get(names[j]), ctx)]);
      }
    }
    for (let round = 0; round < 8; round++) {
      for (const node of internals) {
        let num = 0, den = 0;
        for (const [i, j, m] of pairs) {
          if (m !== node) continue;
          const w = wt(D[i][j]);
          num += w * D[i][j];
          den += 2 * w;
        }
        if (den > 0) heights.set(node, num / den);
      }
      postorder(tree, n => { if (n.children) heights.set(n, Math.max(heights.get(n), ...n.children.map(c => heights.get(c)))); });
    }
    postorder(tree, n => { if (n.children) for (const c of n.children) c.length = Math.max(0, heights.get(n) - heights.get(c)); });
  }

  // ------------------------------------------------------------- parsimony

  function stateMask(v) {
    let m = 0;
    for (let i = 0; i < v.length && i < 31; i++) if (v[i] > 0) m |= (1 << i);
    return m;
  }

  function fitchScore(tree, chars, names) {
    const idx = new Map(names.map((n, i) => [n, i]));
    const leafOf = new Map(leaves(tree).map(l => [l, idx.get(l.name)]));
    let total = 0;
    for (const col of chars) {
      let steps = 0;
      const masks = new Map();
      postorder(tree, n => {
        if (!n.children) { masks.set(n, stateMask(col[leafOf.get(n)])); return; }
        let m = masks.get(n.children[0]);
        for (let k = 1; k < n.children.length; k++) {
          const other = masks.get(n.children[k]);
          const both = m & other;
          if (both) m = both;
          else { m |= other; steps++; }
        }
        masks.set(n, m);
      });
      total += steps;
    }
    return total;
  }

  function sankoffScore(tree, chars, names, cost) {
    const idx = new Map(names.map((n, i) => [n, i]));
    const nStates = cost.length;
    let total = 0;
    for (const col of chars) {
      const dp = new Map();
      postorder(tree, node => {
        if (!node.children) {
          const v = col[idx.get(node.name)];
          const row = new Float64Array(nStates);
          let any = false;
          for (let s = 0; s < nStates; s++) { if (v[s] > 0) { row[s] = 0; any = true; } else row[s] = INF; }
          if (!any) row.fill(0);
          dp.set(node, row);
          return;
        }
        const acc = new Float64Array(nStates);
        for (let s = 0; s < nStates; s++) {
          let sum = 0;
          for (const child of node.children) {
            const cr = dp.get(child);
            let best = INF;
            for (let t = 0; t < nStates; t++) {
              if (cr[t] >= INF) continue;
              const v = cr[t] + (cost[s][t] >= INF ? INF : cost[s][t]);
              if (v < best) best = v;
            }
            sum += best;
          }
          acc[s] = sum >= INF ? INF : sum;
        }
        dp.set(node, acc);
      });
      const r = dp.get(tree);
      const m = Math.min(...r);
      total += isFinite(m) ? m : 0;
    }
    return total;
  }

  function dolloScore(tree, chars, names) {
    const idx = new Map(names.map((n, i) => [n, i]));
    let total = 0;
    for (const col of chars) {
      const dp = new Map();
      postorder(tree, node => {
        if (!node.children) {
          const v = col[idx.get(node.name)];
          const f = [[INF, INF], [INF, INF]];
          const has0 = v[0] > 0;
          let has1 = false;
          for (let st = 1; st < v.length; st++) if (v[st] > 0) { has1 = true; break; }
          if (has1 && !has0) { f[0][1] = 0; f[1][1] = 0; }
          else if (has0 && !has1) { f[0][0] = 0; f[1][0] = 0; }
          else { f[0][0] = 0; f[0][1] = 0; f[1][0] = 0; f[1][1] = 0; }
          dp.set(node, f);
          return;
        }
        const f = [[0, 0], [0, 0]];
        for (let g = 0; g < 2; g++) {
          for (let s = 0; s < 2; s++) {
            let sum = 0;
            for (const child of node.children) {
              const cf = dp.get(child);
              let best = INF;
              for (let c = 0; c < 2; c++) {
                const gain = s === 0 && c === 1 ? 1 : 0;
                if (gain > g) continue;
                const step = s === c ? 0 : 1;
                for (let gc = 0; gc + gain <= g; gc++) {
                  if (cf[gc][c] >= INF) continue;
                  const v = cf[gc][c] + step;
                  if (v < best) best = v;
                }
              }
              sum += best;
            }
            f[g][s] = sum >= INF ? INF : sum;
          }
        }
        dp.set(node, f);
      });
      const r = dp.get(tree);
      const m = Math.min(r[0][0], r[0][1], r[1][0], r[1][1]);
      total += isFinite(m) ? m : 0;
    }
    return total;
  }

  function fitchCost(n) {
    return Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? 0 : 1)));
  }
  function wagnerCost(n) {
    return Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => Math.abs(i - j)));
  }
  function caminSokalCost(n) {
    return Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (j >= i ? j - i : INF)));
  }

  // ------------------------------------------------------------ tree search

  function nniImprove(tree, scoreFn) {
    assignIds(tree);
    let current = scoreFn(tree);
    let anyImprovement = false;
    for (let pass = 0; pass < 30; pass++) {
      const internals = [];
      postorder(tree, n => { if (n.children && n.children.length === 2) internals.push(n); });
      let changed = false;
      for (const node of internals) {
        const candidateNode = findById(tree, node.id);
        if (!candidateNode || !candidateNode.children) continue;
        for (let which = 0; which < 2; which++) {
          const A = candidateNode.children[which];
          if (!A || !A.children || A.children.length !== 2) continue;
          for (let sub = 0; sub < 2; sub++) {
            const trial = cloneTree(tree);
            const tNode = findById(trial, node.id);
            const tA = tNode.children[which];
            const sibling = tNode.children[1 - which];
            const moved = tA.children[sub];
            tA.children[sub] = sibling;
            tNode.children[1 - which] = moved;
            const s = scoreFn(trial);
            if (s > current + 1e-9) {
              tree.children = trial.children;
              tree.name = trial.name;
              current = s;
              changed = true;
              anyImprovement = true;
              break;
            }
          }
          if (changed) break;
        }
        if (changed) break;
      }
      if (!changed) break;
    }
    return anyImprovement;
  }

  /** Stepwise addition driven by a generic tree score (higher is better). */
  function stepwiseSearch(names, scoreOfTree) {
    if (!names.length) return null;
    let tree = { name: null, children: [{ name: names[0], length: 0 }, { name: names[0], length: 0 }] };
    tree = { name: null, children: [{ name: names[0], length: 0 }, { name: names[1], length: 0 }] };
    assignIds(tree);
    for (let k = 2; k < names.length; k++) {
      const taxon = names[k];
      let best = null, bestScore = -Infinity;
      const attachPoints = [];
      (function walk(n, parent) {
        if (parent) attachPoints.push(n.id);
        if (n.children) n.children.forEach(c => walk(c, n));
      })(tree, null);
      for (const id of attachPoints) {
        const trial = cloneTree(tree);
        let placed = false;
        (function walk(n) {
          if (!n.children || placed) return;
          n.children = n.children.map(c => {
            if (c.id === id) {
              placed = true;
              return { name: null, children: [c, { name: taxon, length: 0 }], length: 0 };
            }
            return c;
          });
          n.children.forEach(walk);
        })(trial);
        if (!placed) continue;
        assignIds(trial);
        const s = scoreOfTree(trial);
        if (s > bestScore) { bestScore = s; best = trial; }
      }
      trial_at_root: {
        const trial = { name: null, children: [cloneTree(tree), { name: taxon, length: 0 }], length: null };
        assignIds(trial);
        const s = scoreOfTree(trial);
        if (s > bestScore) { bestScore = s; best = trial; }
      }
      if (!best) break;
      tree = best;
      assignIds(tree);
    }
    return tree;
  }

  // ---------------------------------------------------------- likelihood

  function modelContext(opts) {
    const kind = opts.kind || 'dna';
    const model = PModels.buildQ({
      kind, model: opts.model, freqs: opts.freqs,
      kappa: opts.kappa, kappa2: opts.kappa2, rates: opts.rates
    });
    const rates = opts.alpha && opts.categories > 1 ? PModels.gammaRates(opts.alpha, opts.categories) : [1];
    return { model, rates };
  }

  function tipVector(v, n) {
    const out = new Float64Array(n);
    let any = false;
    for (let i = 0; i < n; i++) if (v[i] > 0) { out[i] = 1; any = true; }
    if (!any) out.fill(1);
    return out;
  }

  function logLikelihood(tree, chars, names, opts, ctx) {
    const c = ctx || modelContext(opts);
    const n = c.model.n;
    const nSites = chars.length;
    const taxonIndex = new Map(names.map((nm, i) => [nm, i]));
    const cache = new Map();
    postorder(tree, node => {
      if (!node.children) {
        const t = taxonIndex.get(node.name);
        const arr = new Array(nSites);
        for (let s = 0; s < nSites; s++) arr[s] = tipVector(chars[s][t], n);
        cache.set(node, arr);
        return;
      }
      const acc = new Array(nSites);
      for (let s = 0; s < nSites; s++) acc[s] = new Float64Array(n).fill(1);
      for (const child of node.children) {
        const cl = cache.get(child);
        const len = Math.max(child.length || 0, 1e-9);
        const combined = new Array(nSites);
        for (let s = 0; s < nSites; s++) combined[s] = new Float64Array(n);
        for (const rate of c.rates) {
          const P = PModels.transitionMatrix(c.model.Q, len * rate, n);
          for (let s = 0; s < nSites; s++) {
            const cv = cl[s], out = combined[s];
            for (let a = 0; a < n; a++) {
              const row = P[a];
              let sum = 0;
              for (let b = 0; b < n; b++) { const v = cv[b]; if (v > 0) sum += row[b] * v; }
              out[a] += sum;
            }
          }
        }
        const inv = 1 / c.rates.length;
        for (let s = 0; s < nSites; s++) {
          const o = combined[s], a2 = acc[s];
          for (let a = 0; a < n; a++) a2[a] *= o[a] * inv;
        }
      }
      cache.set(node, acc);
    });
    const top = cache.get(tree);
    let ll = 0;
    for (let s = 0; s < nSites; s++) {
      let L = 0;
      for (let a = 0; a < n; a++) L += c.model.pi[a] * top[s][a];
      ll += L > 0 ? Math.log(L) : -1e3;
    }
    return ll;
  }

  function optimizeBranchLengths(tree, chars, names, opts, sweeps = 3) {
    const ctx = modelContext(opts);
    const parent = parentMap(tree);
    const objective = () => logLikelihood(tree, chars, names, opts, ctx);
    for (let round = 0; round < sweeps; round++) {
      for (const [child] of parent) {
        const original = child.length;
        let bestT = original && original > 0 ? original : 0.05;
        child.length = bestT;
        let bestV = objective();
        for (const t of [1e-4, 1e-3, 0.01, 0.05, 0.1, 0.2, 0.4, 0.8, 1.5, 3]) {
          child.length = t;
          const v = objective();
          if (v > bestV) { bestV = v; bestT = t; }
        }
        const lo = Math.max(1e-6, bestT / 4), hi = Math.min(6, bestT * 4);
        const fit = PModels.goldenMax(t => { child.length = t; return objective(); }, lo, hi, 18);
        if (fit.fx > bestV) { bestV = fit.fx; bestT = fit.x; }
        child.length = Math.max(1e-6, bestT);
      }
    }
    return objective();
  }

  function pDistances(chars, names, opts) {
    const n = names.length;
    const D = names.map(() => new Array(n).fill(0));
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        let inform = 0, diff = 0;
        for (const col of chars) {
          const a = col[i], b = col[j];
          let dot = 0, sa = 0, sb = 0;
          for (let s = 0; s < a.length; s++) { dot += a[s] * b[s]; sa += a[s]; sb += b[s]; }
          if (!sa || !sb) continue;
          inform++;
          if (dot === 0) diff++;
        }
        const p = inform ? diff / inform : 0;
        const states = chars.length ? chars[0][0].length : 4;
        const b0 = (states - 1) / states;
        const d = p >= b0 ? 2 : -b0 * Math.log(1 - p / b0);
        D[i][j] = D[j][i] = Math.max(1e-6, d);
      }
    }
    return D;
  }

  function mlTree(chars, names, opts = {}) {
    const D = pDistances(chars, names, opts);
    let tree = neighborJoining(D, names);
    assignIds(tree);
    for (const [child] of parentMap(tree)) if (!child.length) child.length = 0.05;
    optimizeBranchLengths(tree, chars, names, opts, 2);
    let current = logLikelihood(tree, chars, names, opts);
    for (let pass = 0; pass < 5; pass++) {
      const before = current;
      nniImprove(tree, t => logLikelihood(t, chars, names, opts));
      optimizeBranchLengths(tree, chars, names, opts, 2);
      current = logLikelihood(tree, chars, names, opts);
      if (current - before < 1e-4) break;
    }
    return { tree, logLikelihood: current, sites: chars.length };
  }

  /** Molecular-clock ML: UPGMA start, then coordinate ascent over node heights. */
  function clockTree(chars, names, opts = {}) {
    const D = pDistances(chars, names, opts);
    const tree = upgma(D, names);
    assignIds(tree);
    const ctx = modelContext(opts);
    const objective = () => logLikelihood(tree, chars, names, opts, ctx);
    const heights = new Map();
    postorder(tree, n => heights.set(n, n.children ? Math.max(...n.children.map(c => heights.get(c) + (c.length || 0))) : 0));
    const sync = () => postorder(tree, n => { if (n.children) for (const c of n.children) c.length = Math.max(1e-6, heights.get(n) - heights.get(c)); });
    sync();
    let best = objective();
    const internals = [];
    postorder(tree, n => { if (n.children) internals.push(n); });
    for (let pass = 0; pass < 4; pass++) {
      for (const node of internals) {
        const floor = Math.max(...node.children.map(c => heights.get(c)));
        const f = h => {
          for (const n of internals) if (n === node) heights.set(n, Math.max(h, floor));
          postorder(tree, n => { if (n.children) heights.set(n, Math.max(heights.get(n), ...n.children.map(c => heights.get(c)))); });
          sync();
          return objective();
        };
        const fit = PModels.goldenMax(f, floor + 1e-4, floor + 3, 18);
        if (fit.fx > best) best = fit.fx;
        heights.set(node, fit.x);
      }
      sync();
      const v = objective();
      if (Math.abs(v - best) < 1e-5) { best = v; break; }
      best = v;
    }
    return { tree, logLikelihood: best, sites: chars.length };
  }

  /** Continuous / gene-frequency data: Brownian-motion distances + least squares. */
  function continuousTree(points, names) {
    const n = names.length;
    const dims = points[0].length;
    const D = names.map(() => new Array(n).fill(0));
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        let num = 0, da = 0, db = 0;
        for (let k = 0; k < dims; k++) {
          num += points[i][k] * points[j][k];
          da += points[i][k] ** 2;
          db += points[j][k] ** 2;
        }
        const corr = da > 0 && db > 0 ? num / Math.sqrt(da * db) : 0;
        D[i][j] = D[j][i] = corr > 0.999999 ? 1e-6 : Math.max(1e-6, -Math.log(Math.max(corr, 1e-6)));
      }
    }
    const fit = leastSquaresFit(D, names, { weights: '1/d2' });
    return { tree: fit.tree, distances: D, dims, ss: fit.ss };
  }

  return {
    upgma, neighborJoining, leastSquaresFit, weightedSS, fitchScore, sankoffScore,
    dolloScore, fitchCost, wagnerCost, caminSokalCost, nniImprove, stepwiseSearch,
    assignIds, findById, logLikelihood, optimizeBranchLengths, mlTree, clockTree,
    continuousTree, pDistances, pathLength, treeIndices, mrca
  };
});
