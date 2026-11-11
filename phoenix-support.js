/* phoenix-support.js — bootstrap resampling, consensus trees, tree comparison,
   invariant site counting, ordination, character compatibility and drawing. */
(function (root, factory) {
  const api = factory(typeof require === 'function' ? require('./phoenix-core.js') : root.PCore);
  if (typeof module === 'object' && module.exports) module.exports = api;
  else root.PSupport = api;
})(typeof globalThis !== 'undefined' ? globalThis : this, function (PCore) {
  'use strict';

  const { leaves, cloneTree, postorder, bipartitions, robinsonFoulds, branchScore, rng } = PCore;

  // -------------------------------------------------------------- bootstrap

  /** Resample the characters (sites) with replacement. */
  function bootstrapReplicates(chars, replicates, seed) {
    const rand = rng(seed === undefined ? 424242 : seed);
    const out = [];
    for (let r = 0; r < replicates; r++) {
      const sample = [];
      for (let k = 0; k < chars.length; k++) sample.push(chars[Math.floor(rand() * chars.length)]);
      out.push(sample);
    }
    return out;
  }

  function jackknifeReplicates(chars, replicates, seed) {
    if (!chars.length) return [];
    const rand = rng(seed === undefined ? 99 : seed);
    const out = [];
    const drop = Math.max(1, Math.round(chars.length * 0.2));
    for (let r = 0; r < replicates; r++) {
      const keep = chars.map((_, i) => i);
      for (let d = 0; d < drop && keep.length > 1; d++) keep.splice(Math.floor(rand() * keep.length), 1);
      out.push(keep.map(i => chars[i]));
    }
    return out;
  }

  // --------------------------------------------------------------- consensus

  function cladeKey(names) {
    return [...names].sort().join('\u0001');
  }

  /** Splits of a tree, each represented by its smaller side (ties broken lexically). */
  function orientedSplits(tree) {
    const all = leaves(tree).map(l => l.name);
    const parent = PCore.parentMap(tree);
    const splits = [];
    const seen = new Set();
    for (const [child] of parent) {
      const below = leaves(child).map(l => l.name);
      if (below.length < 2 || below.length > all.length - 2) continue;
      const complement = all.filter(n => !below.includes(n));
      const inSet = new Set(below);
      const smaller = below.length < complement.length ? below
        : complement.length < below.length ? complement
          : (cladeKey(below) < cladeKey(complement) ? below : complement);
      const key = cladeKey(smaller);
      if (seen.has(key)) continue;
      seen.add(key);
      void inSet;
      splits.push({ key, names: smaller, length: child.length || 0 });
    }
    return splits;
  }

  function areCompatible(a, b) {
    const A = new Set(a.names), B = new Set(b.names);
    let inter = 0;
    for (const x of A) if (B.has(x)) inter++;
    if (inter === 0 || inter === A.size || inter === B.size) return true;
    return false;
  }

  /**
   * Consensus tree from a collection of trees.
   * threshold 0.5 = majority rule, 1.0 = strict consensus.
   */
  function consensus(trees, opts = {}) {
    if (!trees.length) throw new Error('No trees were supplied for a consensus.');
    const threshold = opts.threshold === undefined ? 0.5 : opts.threshold;
    const allTaxa = [...new Set([].concat(...trees.map(t => leaves(t).map(l => l.name))))];
    const counts = new Map();
    const kept = new Map();
    for (const tree of trees) {
      const seen = new Set();
      for (const sp of orientedSplits(tree)) {
        if (seen.has(sp.key)) continue;
        seen.add(sp.key);
        counts.set(sp.key, (counts.get(sp.key) || 0) + 1);
        if (!kept.has(sp.key)) kept.set(sp.key, sp);
      }
    }
    const minimum = threshold > 0 ? Math.ceil(threshold * trees.length - 1e-9) : 1;
    let clades = [...kept.values()]
      .filter(sp => (counts.get(sp.key) || 0) >= minimum)
      .map(sp => ({ names: sp.names, key: sp.key, support: 100 * (counts.get(sp.key) || 0) / trees.length }))
      .sort((a, b) => b.names.length - a.names.length || b.support - a.support);

    // Greedy: keep clades compatible with everything already accepted.
    const accepted = [];
    for (const clade of clades) {
      if (accepted.every(a => areCompatible(a, clade))) accepted.push(clade);
    }
    clades = accepted.sort((a, b) => a.names.length - b.names.length || b.support - a.support);

    const universe = new Set(allTaxa);
    function build(subset, cladeInfo) {
      const key = cladeKey(subset);
      const node = { name: null, children: null, length: 0 };
      if (cladeInfo) node.support = cladeInfo.support;
      const contained = clades.filter(c => c.names.length < subset.size && c.names.length > 1 && c.names.every(n => subset.has(n)));
      const top = contained.filter(c => !contained.some(o => o !== c && o.names.length > c.names.length && c.names.every(n => o.names.includes(n))));
      const covered = new Set();
      const kids = [];
      for (const c of top) {
        const child = build(new Set(c.names), c);
        kids.push(child);
        c.names.forEach(n => covered.add(n));
      }
      for (const n of subset) {
        if (covered.has(n)) continue;
        kids.push({ name: n, children: null, length: 0 });
      }
      void key;
      node.children = kids.length === 1 ? kids : kids;
      if (node.children.length === 1) return node.children[0];
      return node;
    }
    const tree = build(universe, null);
    return { tree, splits: clades, total: trees.length, threshold, allTaxa };
  }

  // ------------------------------------------------------- tree comparison

  function compareTrees(a, b) {
    const rf = robinsonFoulds(a, b);
    const bl = branchScore(a, b);
    return {
      rf: rf.rf,
      maxRf: rf.a + rf.b,
      shared: rf.shared,
      splitsA: rf.a,
      splitsB: rf.b,
      branchScore: bl
    };
  }

  function treeLength(tree) {
    let total = 0;
    PCore.preorder(tree, n => { if (n !== tree) total += n.length || 0; });
    return total;
  }

  // -------------------------------------------------------------- retree ops

  function graphOf(tree) {
    const list = [];
    (function walk(n) {
      const rec = { node: n, id: list.length, edges: [] };
      list.push(rec);
      n.__rec = rec;
      if (n.children) {
        n.children.forEach(c => {
          walk(c);
          const childRec = c.__rec;
          rec.edges.push({ to: childRec, length: c.length || 0 });
          childRec.edges.push({ to: rec, length: c.length || 0 });
        });
      }
    })(tree);
    return { list, byNode: new Map(list.map(r => [r.node, r])) };
  }

  function reroot(tree, outgroup) {
    const g = graphOf(tree);
    const leaf = leaves(tree).find(l => l.name === outgroup);
    if (!leaf) throw new Error(`Outgroup "${outgroup}" is not in the tree.`);
    const start = leaf.__rec;
    if (!start.edges.length) return cloneTree(tree);
    const edge = start.edges[0];
    const build = (rec, from, length) => {
      const kids = rec.edges.filter(e => e.to !== from);
      if (!kids.length) return { name: rec.node.name, children: null, length };
      const node = { name: null, children: kids.map(e => build(e.to, rec, e.length)), length, support: rec.node.support };
      if (node.children.length === 1) {
        const only = node.children[0];
        only.length += length;
        return only;
      }
      return node;
    };
    const side = build(edge.to, start, edge.length / 2);
    return {
      name: null, length: null,
      children: [
        { name: start.node.name, children: null, length: edge.length / 2 },
        side
      ]
    };
  }

  function resolveRandomly(tree, seed) {
    const copy = cloneTree(tree);
    return PCore.resolvePolytomies(copy, rng(seed === undefined ? 1234 : seed));
  }

  function radicalPair(tree, seed) {
    // Randomly flip one internal branch (a "tree rearrangement" operation).
    const copy = cloneTree(tree);
    const rand = rng(seed === undefined ? 77 : seed);
    const internals = [];
    postorder(copy, n => { if (n.children && n.children.length === 2) internals.push(n); });
    const movable = internals.filter(n => n.children.every(c => c.children && c.children.length === 2));
    if (!movable.length) return copy;
    const node = movable[Math.floor(rand() * movable.length)];
    const which = Math.floor(rand() * 2);
    const sub = Math.floor(rand() * 2);
    const A = node.children[which];
    const B = node.children[1 - which];
    const moved = A.children[sub];
    A.children[sub] = B;
    node.children[1 - which] = moved;
    return copy;
  }

  // ---------------------------------------------------------- DNAINVAR core

  /** Per-site invariant support for the three quartet splits. */
  function invariants(chars, names) {
    if (names.length !== 4) throw new Error('DNAINVAR in Phoenix works on exactly four taxa.');
    const counts = { '12|34': 0, '13|24': 0, '14|23': 0 };
    const pairs = { '12|34': [[0, 1], [2, 3]], '13|24': [[0, 2], [1, 3]], '14|23': [[0, 3], [1, 2]] };
    let informative = 0, uninformative = 0;
    const patternTotals = new Map();
    for (const col of chars) {
      const code = col.map(v => {
        let best = -1, bestCount = 0;
        for (let s = 0; s < v.length; s++) if (v[s] > bestCount) { bestCount = v[s]; best = s; }
        let alternatives = 0;
        for (let s = 0; s < v.length; s++) if (v[s] > 0) alternatives++;
        return alternatives === 1 ? best : -1;
      });
      if (code.some(c => c < 0)) { uninformative++; continue; }
      const key = code.join('');
      patternTotals.set(key, (patternTotals.get(key) || 0) + 1);
      let supported = 0;
      for (const [split, [[a, b], [c, d]]] of Object.entries(pairs)) {
        if (code[a] === code[b] && code[c] === code[d] && code[a] !== code[c]) { counts[split]++; supported++; }
      }
      if (supported) informative++;
      else uninformative++;
    }
    const best = Object.entries(counts).sort((x, y) => y[1] - x[1])[0];
    const total = Object.values(counts).reduce((a, b) => a + b, 0);
    return {
      counts, informative, uninformative, patterns: patternTotals,
      bestSplit: best && best[1] > 0 ? best[0] : null,
      bestCount: best ? best[1] : 0,
      ratio: total ? best[1] / total : 0,
      total
    };
  }

  // -------------------------------------------------------------- ordination

  /** Principal coordinates (Gower) of a distance matrix. */
  function pcoa(D, k = 3) {
    const n = D.length;
    const d2 = D.map(row => row.map(v => v * v));
    const rowMean = d2.map(row => row.reduce((a, b) => a + b, 0) / n);
    const grand = rowMean.reduce((a, b) => a + b, 0) / n;
    const B = [];
    for (let i = 0; i < n; i++) {
      B.push([]);
      for (let j = 0; j < n; j++) B[i].push(-0.5 * (d2[i][j] - rowMean[i] - rowMean[j] + grand));
    }
    const { values, vectors } = jacobiEigen(B);
    const order = values.map((v, i) => i).sort((a, b) => values[b] - values[a]);
    const axes = [];
    for (let a = 0; a < Math.min(k, n); a++) {
      const idx = order[a];
      const scale = Math.sqrt(Math.max(values[idx], 1e-12));
      axes.push(vectors.map(row => row[idx] * scale));
    }
    const positive = values.filter(v => v > 1e-9);
    const totalVar = positive.reduce((a, b) => a + b, 0);
    const explained = axes.map((_, a) => (positive[order[a]] || 0) / totalVar);
    return { axes, eigenvalues: order.map(i => values[i]), explained };
  }

  /** Jacobi eigenvalue algorithm for symmetric matrices. */
  function jacobiEigen(A) {
    const n = A.length;
    const a = A.map(r => r.slice());
    let v = [];
    for (let i = 0; i < n; i++) { v.push(new Array(n).fill(0)); v[i][i] = 1; }
    for (let sweep = 0; sweep < 100; sweep++) {
      let off = 0;
      for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) off += a[i][j] * a[i][j];
      if (off < 1e-18) break;
      for (let p = 0; p < n - 1; p++) {
        for (let q = p + 1; q < n; q++) {
          if (Math.abs(a[p][q]) < 1e-14) continue;
          const theta = (a[q][q] - a[p][p]) / (2 * a[p][q]);
          const t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
          const c = 1 / Math.sqrt(t * t + 1), s = t * c;
          for (let i = 0; i < n; i++) {
            const aip = a[i][p], aiq = a[i][q];
            a[i][p] = c * aip - s * aiq;
            a[i][q] = s * aip + c * aiq;
          }
          for (let i = 0; i < n; i++) {
            const api = a[p][i], aqi = a[q][i];
            a[p][i] = c * api - s * aqi;
            a[q][i] = s * api + c * aqi;
          }
          for (let i = 0; i < n; i++) {
            const vip = v[i][p], viq = v[i][q];
            v[i][p] = c * vip - s * viq;
            v[i][q] = s * vip + c * viq;
          }
        }
      }
    }
    return { values: a.map((row, i) => row[i]), vectors: v };
  }

  // ------------------------------------------------ character compatibility

  function characterStates(chars) {
    // chars[site][taxon] -> state vector; return per-site map taxon->state code
    return chars.map(col => col.map(v => {
      let best = 0, bestCount = -1, count = 0;
      for (let s = 0; s < v.length; s++) {
        if (v[s] > 0) { count++; if (v[s] > bestCount) { bestCount = v[s]; best = s; } }
      }
      return count === 1 ? best : -1;
    }));
  }

  function compatible(siteA, siteB) {
    const sets = new Map();
    for (let t = 0; t < siteA.length; t++) {
      const a = siteA[t], b = siteB[t];
      if (a < 0 || b < 0) continue;
      if (!sets.has(a)) sets.set(a, new Set());
      sets.get(a).add(b);
    }
    const entries = [...sets.values()];
    for (let i = 0; i < entries.length; i++) {
      for (let j = i + 1; j < entries.length; j++) {
        const X = entries[i], Y = entries[j];
        let inter = 0;
        for (const x of X) if (Y.has(x)) inter++;
        if (inter === 0) continue;
        if (inter === X.size || inter === Y.size) continue;
        return false;
      }
    }
    return true;
  }

  /** Largest set of mutually compatible characters (Bron–Kerbosch). */
  function cliques(chars) {
    const sites = characterStates(chars);
    const n = sites.length;
    const adj = sites.map(() => new Array(n).fill(false));
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        const ok = compatible(sites[i], sites[j]);
        adj[i][j] = ok; adj[j][i] = ok;
      }
    }
    let best = [];
    const expand = (R, P) => {
      if (!P.length) { if (R.length > best.length) best = [...R]; return; }
      if (R.length + P.length <= best.length) return;
      for (let k = 0; k < P.length; k++) {
        const v = P[k];
        const next = P.slice(k + 1).filter(w => adj[v][w]);
        expand([...R, v], next);
      }
    };
    expand([], sites.map((_, i) => i));
    return { indices: best, count: best.length, totalCharacters: n };
  }

  // --------------------------------------------------------------- rendering

  function formatLength(v) {
    if (!isFinite(v)) return '—';
    if (Math.abs(v) >= 0.001) return v.toFixed(4);
    return v.toExponential(2);
  }

  /** Rectangular phylogram/cladogram as SVG. */
  function drawRectangular(tree, opts = {}) {
    const names = leaves(tree).map(l => l.name);
    const width = opts.width || 720;
    const rowHeight = opts.rowHeight || 22;
    const padTop = 24, padBottom = 20, padLeft = 20, labelWidth = opts.labelWidth || 190;
    const height = Math.max(160, names.length * rowHeight + padTop + padBottom);
    const plotWidth = width - padLeft - labelWidth;
    const cladogram = opts.style === 'cladogram';
    const ctx = buildLayout(tree, cladogram);
    const maxX = ctx.maxX || 1;
    const y = new Map();
    names.forEach((n, i) => y.set(n, padTop + i * rowHeight + rowHeight / 2));
    const sx = x => padLeft + (x / maxX) * plotWidth;
    const parts = [];
    function assign(node) {
      if (!node.children) { node.__y = y.get(node.name); node.__x = sx(node.__depth || 0); return; }
      node.children.forEach(assign);
      node.__y = (node.children[0].__y + node.children[node.children.length - 1].__y) / 2;
      node.__x = sx(node.__depth || 0);
    }
    assign(tree);
    const lines = [];
    function walk(node) {
      if (!node.children) return;
      for (const c of node.children) {
        lines.push(`<path d="M${node.__x.toFixed(2)},${node.__y.toFixed(2)} H${c.__x.toFixed(2)} V${c.__y.toFixed(2)}" fill="none" stroke="#7d9c84" stroke-width="1.6" stroke-linecap="round"/>`);
        walk(c);
      }
    }
    walk(tree);
    const labels = names.map(n => {
      const leaf = leaves(tree).find(l => l.name === n);
      const yy = y.get(n);
      const label = n.length > 26 ? n.slice(0, 25) + '…' : n;
      return `<text x="${(sx(leaf.__depth || 0) + 10).toFixed(2)}" y="${(yy + 4).toFixed(2)}" font-family="DM Sans, system-ui, sans-serif" font-size="12" fill="#3f5d47">${escapeXml(label)}</text>`;
    }).join('');
    const ticks = cladogram ? '' : axisTicks(maxX, sx, padTop, height - padBottom);
    parts.push(`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" role="img" aria-label="Phylogenetic tree" style="width:100%;height:auto">`);
    parts.push(ticks);
    parts.push(lines.join(''));
    parts.push(labels);
    parts.push('</svg>');
    return parts.join('');
  }

  function axisTicks(maxX, sx, top, bottom) {
    const out = [`<line x1="${sx(0)}" y1="${top - 8}" x2="${sx(maxX)}" y2="${top - 8}" stroke="#d8e2d5" stroke-width="1"/>`];
    for (let i = 0; i <= 4; i++) {
      const v = (maxX * i) / 4;
      out.push(`<line x1="${sx(v)}" y1="${top - 12}" x2="${sx(v)}" y2="${top - 4}" stroke="#c9d7c5"/>`);
      out.push(`<text x="${sx(v)}" y="${top - 16}" text-anchor="middle" font-family="DM Sans, system-ui, sans-serif" font-size="9" fill="#93a894">${v.toFixed(3)}</text>`);
    }
    void bottom;
    return out.join('');
  }

  function buildLayout(tree, cladogram) {
    let maxX = 0;
    (function walk(node, depth) {
      node.__depth = cladogram ? depth : depth === 0 ? 0 : node.__parentDepth + (node.length || 0);
      maxX = Math.max(maxX, node.__depth || 0);
      if (node.children) node.children.forEach(c => { c.__parentDepth = node.__depth; walk(c, depth + 1); });
    })(tree, 0);
    if (cladogram) maxX = maxX || 1;
    return { maxX };
  }

  /** Radial (unrooted-looking) tree for DRAWTREE. */
  function drawRadial(tree, opts = {}) {
    const names = leaves(tree).map(l => l.name);
    const size = opts.size || 620;
    const cx = size / 2, cy = size / 2;
    const radius = size / 2 - (opts.margin || 90);
    const ctx = buildLayout(tree, opts.style === 'cladogram');
    const maxX = ctx.maxX || 1;
    const angle = new Map();
    names.forEach((n, i) => angle.set(n, (2 * Math.PI * i) / names.length - Math.PI / 2));
    function place(node) {
      if (!node.children) { node.__a = angle.get(node.name); node.__r = (node.__depth || 0) / maxX; return; }
      node.children.forEach(place);
      const a = node.children.map(c => c.__a);
      let sinSum = 0, cosSum = 0;
      for (const x of a) { sinSum += Math.sin(x); cosSum += Math.cos(x); }
      node.__a = Math.atan2(sinSum / a.length, cosSum / a.length);
      node.__r = (node.__depth || 0) / maxX;
    }
    place(tree);
    const P = (a, r) => [cx + Math.cos(a) * r * radius, cy + Math.sin(a) * r * radius];
    const seg = [];
    (function walk(node) {
      if (!node.children) return;
      for (const c of node.children) {
        const [x1, y1] = P(node.__a, node.__r);
        const [x2, y2] = P(c.__a, c.__r);
        seg.push(`<line x1="${x1.toFixed(1)}" y1="${y1.toFixed(1)}" x2="${x2.toFixed(1)}" y2="${y2.toFixed(1)}" stroke="#7d9c84" stroke-width="1.6" stroke-linecap="round"/>`);
        walk(c);
      }
    })(tree);
    const labels = names.map(n => {
      const leaf = leaves(tree).find(l => l.name === n);
      const [x, y] = P(leaf.__a, Math.min(1.06, leaf.__r + 0.04));
      const anchor = Math.cos(leaf.__a) >= 0 ? 'start' : 'end';
      const label = n.length > 20 ? n.slice(0, 19) + '…' : n;
      return `<text x="${x.toFixed(1)}" y="${(y + 4).toFixed(1)}" text-anchor="${anchor}" font-family="DM Sans, system-ui, sans-serif" font-size="11" fill="#3f5d47">${escapeXml(label)}</text>`;
    }).join('');
    return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${size} ${size}" role="img" aria-label="Unrooted phylogram" style="width:100%;height:auto">${seg.join('')}${labels}</svg>`;
  }

  function escapeXml(s) {
    return String(s).replace(/[<>&'"]/g, c => ({ '<': '&lt;', '>': '&gt;', '&': '&amp;', "'": '&apos;', '"': '&quot;' }[c]));
  }

  /** Plain-text tree for the report pane. */
  function asciiTree(tree, opts = {}) {
    const showLengths = opts.lengths !== false;
    const lines = [];
    const label = n => {
      const name = n.name || '';
      if (!showLengths || n.length === undefined || n.length === null) return name;
      return name ? `${name}:${formatLength(n.length)}` : `${formatLength(n.length)}`;
    };
    (function walk(node, prefix, isLast, isRoot) {
      const connector = isRoot ? '' : (isLast ? '└── ' : '├── ');
      const text = label(node);
      if (!(isRoot && !text)) lines.push(prefix + connector + text);
      if (node.children) {
        node.children.forEach((c, i) => {
          const last = i === node.children.length - 1;
          walk(c, prefix + (isRoot ? '' : (isLast ? '    ' : '│   ')), last, false);
        });
      }
    })(tree, '', true, true);
    return lines.join('\n');
  }

  return {
    bootstrapReplicates, jackknifeReplicates, consensus, orientedSplits,
    compareTrees, treeLength, reroot, resolveRandomly, radicalPair,
    invariants, pcoa, jacobiEigen, compatible, cliques,
    drawRectangular, drawRadial, asciiTree, escapeXml, formatLength
  };
});
