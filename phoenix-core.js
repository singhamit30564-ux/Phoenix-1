/* phoenix-core.js — data model, parsing and tree primitives.
   Shared by every PHYLIP-style program in the suite. */
(function (root, factory) {
  const api = factory();
  if (typeof module === 'object' && module.exports) module.exports = api;
  else root.PCore = api;
})(typeof globalThis !== 'undefined' ? globalThis : this, function () {
  'use strict';

  const DNA = ['A', 'C', 'G', 'T'];
  const AA = 'ARNDCQEGHILKMFPSTWYV'.split('');
  const BIN = ['0', '1'];

  const DNA_AMBIG = {
    A: ['A'], C: ['C'], G: ['G'], T: ['T'], U: ['T'],
    R: ['A', 'G'], Y: ['C', 'T'], S: ['C', 'G'], W: ['A', 'T'],
    K: ['G', 'T'], M: ['A', 'C'], B: ['C', 'G', 'T'], D: ['A', 'G', 'T'],
    H: ['A', 'C', 'T'], V: ['A', 'C', 'G']
  };
  const AA_AMBIG = { B: ['D', 'N'], Z: ['E', 'Q'], J: ['I', 'L'] };
  const MISSING_CHARS = 'N?-.X*';

  /** Character state alphabet for a program's input kind. */
  function alphabet(kind) {
    if (kind === 'protein') return { kind, states: AA, ambig: AA_AMBIG, missing: 'X?-.B*ZJ' };
    if (kind === 'sites') return { kind, states: BIN, ambig: {}, missing: '?' };
    if (kind === 'discrete') return { kind, states: '0123456789'.split(''), ambig: {}, missing: '?' };
    return { kind: 'dna', states: DNA, ambig: DNA_AMBIG, missing: 'N?-.X*' };
  }

  function stateVector(ch, alpha) {
    const up = ch.toUpperCase();
    const n = alpha.states.length;
    const out = new Float64Array(n);
    if (alpha.missing.includes(up) && !alpha.ambig[up]) { out.fill(1); return out; }
    const list = alpha.ambig[up] || (alpha.states.includes(up) ? [up] : null);
    if (!list) throw new Error(`"${ch}" is not a valid ${alpha.kind} character.`);
    for (const s of list) {
      const idx = alpha.states.indexOf(s);
      if (idx >= 0) out[idx] = 1;
    }
    return out;
  }

  function isMissing(ch, alpha) {
    const up = ch.toUpperCase();
    if (alpha.ambig[up]) return false;
    return alpha.missing.includes(up) || !alpha.states.includes(up);
  }

  // ---------------------------------------------------------------- parsing

  /** Parse FASTA. Returns [{name, seq}] with whitespace stripped. */
  function parseFasta(text) {
    const records = [];
    let cur = null;
    for (const raw of String(text).split(/\r?\n/)) {
      const line = raw.trim();
      if (!line || line.startsWith(';')) continue;
      if (line[0] === '>') {
        const name = line.slice(1).trim().split(/\s+/)[0];
        if (!name) throw new Error('Every FASTA header needs a name.');
        cur = { name, seq: '' };
        records.push(cur);
      } else {
        if (!cur) throw new Error('Sequences must begin with a FASTA header (>name).');
        cur.seq += line.replace(/\s+/g, '').toUpperCase();
      }
    }
    return records;
  }

  /** Parse sequential or interleaved PHYLIP format. */
  function parsePhylip(text) {
    const lines = String(text).split(/\r?\n/).filter(l => l.trim().length);
    if (!lines.length) throw new Error('PHYLIP input is empty.');
    const head = lines[0].trim().split(/\s+/);
    const n = parseInt(head[0], 10);
    const len = parseInt(head[1], 10);
    if (!n || !len || n < 1 || len < 1) throw new Error('First line must be "<sequences> <sites>".');
    const names = [], seqs = [];
    let i = 1;
    for (let k = 0; k < n; k++) {
      const line = lines[i++];
      if (line === undefined) throw new Error(`Expected ${n} sequences, found ${k}.`);
      const name = line.slice(0, 10).trim() || line.trim().split(/\s+/)[0];
      names.push(name);
      seqs.push(line.slice(10).replace(/\s+/g, '').toUpperCase());
    }
    // Interleaved continuation blocks.
    while (i < lines.length && seqs.some(s => s.length < len)) {
      for (let k = 0; k < n && i < lines.length; k++) {
        const line = lines[i++];
        if (!line) continue;
        seqs[k] += line.replace(/\s+/g, '').toUpperCase();
      }
    }
    const bad = seqs.findIndex(s => s.length !== len);
    if (bad >= 0) throw new Error(`Sequence "${names[bad]}" has ${seqs[bad].length} sites, expected ${len}.`);
    return names.map((name, k) => ({ name, seq: seqs[k] }));
  }

  /** Auto-detect FASTA vs PHYLIP and parse. */
  function parseInput(text) {
    const body = String(text).trim();
    if (!body) throw new Error('No input data provided.');
    return body[0] === '>' ? parseFasta(body) : parsePhylip(body);
  }

  function validateRecords(records, alpha, opts = {}) {
    const max = opts.max || 200;
    if (records.length < (opts.min || 2)) throw new Error(`At least ${opts.min || 2} sequences are required; got ${records.length}.`);
    if (records.length > max) throw new Error(`At most ${max} taxa are supported by this program in Phoenix.`);
    if (new Set(records.map(r => r.name)).size !== records.length) throw new Error('Taxon names must be unique.');
    const len = records[0].seq.length;
    if (!len) throw new Error('Sequences are empty.');
    if (records.some(r => r.seq.length !== len)) throw new Error('All sequences must be aligned to the same length.');
    for (const r of records) {
      for (const ch of r.seq) {
        const up = ch.toUpperCase();
        if (!alpha.states.includes(up) && !alpha.ambig[up] && !alpha.missing.includes(up)) {
          throw new Error(`Character "${ch}" in "${r.name}" is not valid for ${alpha.kind} input.`);
        }
      }
    }
    return { length: len };
  }

  /** Build per-site character vectors: chars[site][taxon] -> state vector. */
  function encode(records, alpha) {
    const len = records[0].seq.length;
    const out = [];
    for (let s = 0; s < len; s++) {
      const col = records.map(r => stateVector(r.seq[s], alpha));
      // A site where nothing is observed carries no information.
      if (col.every(v => v.every(x => x === 0))) continue;
      out.push(col);
    }
    return out;
  }

  /** Observed state index per taxon, -1 when the site is unobservable there. */
  function observedCodes(records, alpha) {
    const len = records[0].seq.length;
    const out = [];
    for (let s = 0; s < len; s++) {
      const col = records.map(r => {
        const ch = r.seq[s].toUpperCase();
        if (isMissing(ch, alpha) || alpha.ambig[ch]) return -1;
        return alpha.states.indexOf(ch);
      });
      out.push({ site: s, codes: col });
    }
    return out;
  }

  // ---------------------------------------------------------------- newick

  function parseNewick(text) {
    const s = String(text).trim();
    let i = 0;
    const ws = () => { while (i < s.length && /\s/.test(s[i])) i++; };
    function subtree() {
      ws();
      const node = { children: null, name: null, length: null, support: null };
      if (s[i] === '(') {
        i++; node.children = [];
        for (;;) {
          node.children.push(subtree());
          ws();
          if (s[i] === ',') { i++; continue; }
          if (s[i] === ')') { i++; break; }
          throw new Error('Malformed Newick: expected "," or ")".');
        }
        ws();
        let label = '';
        while (i < s.length && !',):;'.includes(s[i])) label += s[i++];
        label = label.trim();
        const num = Number(label);
        if (label && !Number.isNaN(num) && /^[0-9.]+$/.test(label)) node.support = num;
        else if (label) node.name = label.replace(/^['"]|['"]$/g, '');
      } else {
        let label = '';
        while (i < s.length && !',():;'.includes(s[i])) label += s[i++];
        label = label.trim().replace(/^['"]|['"]$/g, '');
        if (!label) throw new Error('Malformed Newick: missing taxon name.');
        node.name = label;
      }
      ws();
      if (s[i] === ':') {
        i++; let num = '';
        while (i < s.length && !',);'.includes(s[i])) num += s[i++];
        const v = parseFloat(num);
        if (Number.isNaN(v)) throw new Error(`Malformed Newick: bad branch length "${num}".`);
        node.length = v;
      }
      return node;
    }
    const tree = subtree();
    return tree;
  }

  function quoteName(name) {
    const safe = /^[A-Za-z0-9_.-]+$/.test(name);
    return safe ? name : "'" + String(name).replace(/'/g, "''") + "'";
  }

  function toNewick(node, opts = {}) {
    const decimals = opts.decimals === undefined ? 6 : opts.decimals;
    const withSupport = opts.support !== false;
    function walk(n) {
      let label;
      if (n.children && n.children.length) {
        label = '(' + n.children.map(walk).join(',') + ')';
        if (withSupport && n.support !== undefined && n.support !== null) label += n.support.toFixed(opts.supportDecimals === undefined ? 1 : opts.supportDecimals);
        if (n.name) label += quoteName(n.name);
      } else {
        label = quoteName(n.name === null || n.name === undefined ? '' : n.name);
      }
      if (n.length !== undefined && n.length !== null && !Number.isNaN(n.length)) {
        const v = Math.abs(n.length) < 1e-12 ? 0 : n.length;
        label += ':' + v.toFixed(decimals);
      }
      return label;
    }
    return walk(node) + ';';
  }

  // ---------------------------------------------------------------- trees

  function leaves(node, out = []) {
    if (!node.children || !node.children.length) out.push(node);
    else node.children.forEach(c => leaves(c, out));
    return out;
  }

  function cloneTree(node) {
    return {
      id: node.id, name: node.name, length: node.length, support: node.support,
      height: node.height, size: node.size,
      children: node.children ? node.children.map(cloneTree) : null
    };
  }

  function numberNodes(node) {
    let id = 0;
    (function walk(n) { n.id = id++; if (n.children) n.children.forEach(walk); })(node);
    return node;
  }

  function postorder(node, fn) {
    if (node.children) node.children.forEach(c => postorder(c, fn));
    fn(node);
  }

  function preorder(node, fn) {
    fn(node);
    if (node.children) node.children.forEach(c => preorder(c, fn));
  }

  function parentMap(root) {
    const map = new Map();
    preorder(root, n => { if (n.children) n.children.forEach(c => map.set(c, n)); });
    return map;
  }

  /** Canonical bipartition key (leaf-name set side with fewest members). */
  function splitKey(names) {
    const a = [...names].sort();
    const b = names._all.filter(x => !names.has(x)).sort();
    const side = a.length < b.length || (a.length === b.length && a[0] < b[0]) ? a : b;
    return side.join('\u0001');
  }

  /** Set of internal bipartitions of an unrooted/rooted tree. */
  function bipartitions(root) {
    const all = leaves(root).map(l => l.name);
    const parent = parentMap(root);
    const keys = new Set();
    for (const [child, par] of parent) {
      const below = new Set(leaves(child).map(l => l.name));
      if (below.size < 2 || below.size > all.length - 2) continue;
      const a = [...below].sort().join('\u0001');
      const b = all.filter(n => !below.has(n)).sort().join('\u0001');
      keys.add(a < b ? a : b);
    }
    return keys;
  }

  /** Robinson–Foulds symmetric difference between two trees. */
  function robinsonFoulds(a, b) {
    const A = bipartitions(a), B = bipartitions(b);
    let shared = 0;
    for (const k of A) if (B.has(k)) shared++;
    return { rf: (A.size - shared) + (B.size - shared), shared, a: A.size, b: B.size };
  }

  function branchScore(a, b) {
    const collect = root => {
      const parent = parentMap(root);
      const m = new Map();
      for (const [child] of parent) {
        const below = new Set(leaves(child).map(l => l.name));
        const all = leaves(root).map(l => l.name);
        if (below.size < 1 || below.size > all.length - 1) continue;
        const a = [...below].sort().join('\u0001');
        const b = all.filter(n => !below.has(n)).sort().join('\u0001');
        const key = a < b ? a : b;
        m.set(key, (m.get(key) || 0) + Math.abs(child.length || 0));
      }
      return m;
    };
    const A = collect(a), B = collect(b);
    let sum = 0;
    const keys = new Set([...A.keys(), ...B.keys()]);
    for (const k of keys) {
      const d = (A.get(k) || 0) - (B.get(k) || 0);
      sum += d * d;
    }
    return Math.sqrt(sum);
  }

  function resolvePolytomies(node, rng) {
    if (!node.children) return node;
    node.children.forEach(c => resolvePolytomies(c, rng));
    while (node.children.length > 2) {
      const j = 1 + Math.floor(rng() * (node.children.length - 1));
      const picked = [node.children[0], node.children[j]];
      node.children = node.children.filter(c => !picked.includes(c));
      node.children.push({ children: picked, name: null, length: null });
    }
    return node;
  }

  function distanceMatrix(names, fn) {
    const m = names.map(() => new Array(names.length).fill(0));
    for (let i = 0; i < names.length; i++) {
      for (let j = i + 1; j < names.length; j++) {
        const d = fn(i, j);
        m[i][j] = d; m[j][i] = d;
      }
    }
    return m;
  }

  function formatMatrix(names, m) {
    const head = '      ' + names.map((_, j) => String(j + 1).padStart(9)).join('');
    const rows = m.map((row, i) => String(i + 1).padStart(5) + row.map(v => v.toFixed(5).padStart(9)).join(''));
    return [head, ...rows].join('\n');
  }

  function rng(seed) {
    let s = (seed === undefined ? 12345 : seed) >>> 0;
    return function () {
      s ^= s << 13; s >>>= 0;
      s ^= s >> 17;
      s ^= s << 5; s >>>= 0;
      return s / 4294967296;
    };
  }

  return {
    DNA, AA, BIN, alphabet, stateVector, isMissing, parseFasta, parsePhylip, parseInput,
    validateRecords, encode, observedCodes, parseNewick, toNewick, quoteName,
    leaves, cloneTree, numberNodes, postorder, preorder, parentMap, bipartitions,
    robinsonFoulds, branchScore, resolvePolytomies, distanceMatrix, formatMatrix, rng
  };
});
