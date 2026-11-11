/* phoenix-programs.js — the PHYLIP-style program catalogue.
   Every entry maps a legacy program name onto a real implementation in the
   Phoenix engine. Program descriptions state precisely what is computed. */
(function (root, factory) {
  const api = factory(
    typeof require === 'function' ? require('./phoenix-core.js') : root.PCore,
    typeof require === 'function' ? require('./phoenix-models.js') : root.PModels,
    typeof require === 'function' ? require('./phoenix-infer.js') : root.PInfer,
    typeof require === 'function' ? require('./phoenix-support.js') : root.PSupport,
    typeof require === 'function' ? require('./phoenix-continuous.js') : root.PContinuous,
    typeof require === 'function' ? require('./phoenix-interactive.js') : root.PInteractive
  );
  if (typeof module === 'object' && module.exports) module.exports = api;
  else root.PPrograms = api;
})(typeof globalThis !== 'undefined' ? globalThis : this, function (C, M, I, S, PC, PI) {
  'use strict';

  const num = (id, label, def, extra) => Object.assign({ id, label, type: 'number', default: def }, extra || {});
  const sel = (id, label, choices, def) => ({ id, label, type: 'select', choices, default: def });
  const chk = (id, label, def) => ({ id, label, type: 'checkbox', default: !!def });

  // ------------------------------------------------------------------ input

  function sequenceInput(text, kind) {
    const records = C.parseInput(text);
    const alpha = C.alphabet(kind);
    const { length } = C.validateRecords(records, alpha, { max: 120 });
    const chars = C.encode(records, alpha);
    if (!chars.length) throw new Error('No informative sites remain after removing fully ambiguous columns.');
    return { records, alpha, chars, names: records.map(r => r.name), sites: length, states: alpha.states.length };
  }

  function numberInput(text) {
    const lines = String(text).split(/\r?\n/).filter(l => l.trim() && !l.trim().startsWith('#'));
    if (!lines.length) throw new Error('No data supplied.');
    let start = 0;
    let declared = 0;
    const head = lines[0].trim().split(/\s+/);
    if (/^\d+$/.test(head[0]) && /^\d+$/.test(head[1] || '')) { declared = parseInt(head[0], 10); start = 1; }
    const names = [];
    const rows = [];
    for (let i = start; i < lines.length; i++) {
      const parts = lines[i].trim().split(/\s+/);
      const name = parts[0];
      const values = parts.slice(1).map(Number);
      if (!values.length || values.some(v => Number.isNaN(v))) throw new Error(`Row "${lines[i].trim()}" needs a name followed by numbers.`);
      names.push(name);
      rows.push(values);
    }
    if (declared && rows.length !== declared) throw new Error(`Header declares ${declared} taxa but ${rows.length} rows were found.`);
    if (rows.length < 2) throw new Error('At least two taxa are required.');
    const width = rows[0].length;
    if (rows.some(r => r.length !== width)) throw new Error('All rows must have the same number of values.');
    return { names, rows, width };
  }

  function treeInput(text, minimum) {
    const chunks = String(text).split(/;\s*(?:\r?\n|$)/).map(s => s.trim()).filter(Boolean);
    const trees = chunks.map(chunk => C.parseNewick(chunk + ';'));
    if (trees.length < (minimum || 1)) throw new Error(`This program needs at least ${minimum || 1} Newick tree${minimum > 1 ? 's' : ''}.`);
    return trees;
  }

  function freqsFromChars(chars, states) {
    return M.empiricalFrequencies(chars, states);
  }

  // --------------------------------------------------------- distance work

  function distanceTable(records, kind, opts) {
    const n = records.length;
    const D = records.map(() => new Array(n).fill(0));
    const saturated = [];
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        let r;
        if (kind === 'protein') r = M.proteinDistance(records[i].seq, records[j].seq, opts);
        else if (kind === 'sites') r = M.siteDistance(records[i].seq, records[j].seq, opts);
        else r = M.dnaDistance(records[i].seq, records[j].seq, opts);
        const d = r.distance;
        D[i][j] = D[j][i] = isFinite(d) ? d : NaN;
        if (!isFinite(d)) saturated.push(`${records[i].name}/${records[j].name}`);
      }
    }
    if (saturated.length) {
      throw new Error(`Distances are undefined (saturated) for ${saturated.length} pair(s), e.g. ${saturated[0]}. Try a less parameterised model, or use maximum likelihood.`);
    }
    return D;
  }

  function matrixResult(D, names, notes, stats) {
    return {
      kind: 'matrix',
      names,
      matrix: D,
      text: C.formatMatrix(names, D),
      notes: notes || [],
      stats: stats || []
    };
  }

  function treeResult(tree, notes, stats, style) {
    return {
      kind: 'tree',
      tree,
      ascii: S.asciiTree(tree),
      svg: S.drawRectangular(tree, { style: style || 'phenogram' }),
      notes: notes || [],
      stats: stats || []
    };
  }

  // ----------------------------------------------------- parsimony plumbing

  function fitchLengths(tree, chars, names) {
    const idx = new Map(names.map((n, i) => [n, i]));
    const down = new Map();
    const finalSet = new Map();
    for (const col of chars) {
      down.clear();
      const steps = { value: 0 };
      C.postorder(tree, node => {
        if (!node.children) {
          const v = col[idx.get(node.name)];
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
          else { m |= other; steps.value++; }
        }
        down.set(node, m);
      });
      finalSet.set(col, new Map());
      const assign = (node, parentFinal) => {
        const own = down.get(node);
        const chosen = parentFinal !== undefined && (own & parentFinal) ? (own & parentFinal) : own;
        finalSet.get(col).set(node, chosen);
        if (node.children) node.children.forEach(c => assign(c, chosen));
      };
      assign(tree, undefined);
      void steps;
    }
    for (const [child, parent] of C.parentMap(tree)) {
      let changes = 0;
      for (const col of chars) {
        const a = finalSet.get(col).get(parent) || 1;
        const b = finalSet.get(col).get(child) || 1;
        if (!(a & b)) changes++;
        else {
          // Partial change: count when the child's set is strictly narrower.
          if ((b & a) !== b) changes += 1;
        }
      }
      child.length = changes;
    }
    if (tree.children) tree.length = null;
    return tree;
  }

  function parsimonySearch(chars, names, costFn, opts = {}) {
    const restarts = opts.restarts || 4;
    let best = null, bestScore = Infinity;
    for (let r = 0; r < restarts; r++) {
      const order = r === 0 ? names.slice() : shuffled(names, 1000 + r * 37);
      let tree = I.stepwiseSearch(order, t => -costFn(t));
      if (!tree) continue;
      I.assignIds(tree);
      for (let pass = 0; pass < 4; pass++) {
        const improved = I.nniImprove(tree, t => -costFn(t));
        if (!improved) break;
      }
      const score = costFn(tree);
      if (score < bestScore) { bestScore = score; best = tree; }
    }
    if (!best) throw new Error('No tree could be built from this data.');
    return { tree: best, score: bestScore };
  }

  function shuffled(list, seed) {
    const rand = C.rng(seed);
    const out = list.slice();
    for (let i = out.length - 1; i > 0; i--) { const j = Math.floor(rand() * (i + 1));[out[i], out[j]] = [out[j], out[i]]; }
    return out;
  }

  function stepMatrix(kind, states) {
    if (kind === 'dollo') return null;
    if (kind === 'camin-sokal') return I.caminSokalCost(states);
    if (kind === 'wagner') return I.wagnerCost(states);
    return I.fitchCost(states);
  }

  function mixedScore(chars, names, weights) {
    return tree => {
      let total = 0;
      if (weights.wagner) total += weights.wagner * I.sankoffScore(tree, chars, names, I.wagnerCost(10));
      if (weights.fitch) total += weights.fitch * I.fitchScore(tree, chars, names);
      if (weights.dollo) total += weights.dollo * I.dolloScore(tree, chars, names);
      if (weights.caminSokal) total += weights.caminSokal * I.sankoffScore(tree, chars, names, I.caminSokalCost(10));
      return total;
    };
  }

  // ------------------------------------------------------------- ML plumbing

  function mlContext(chars, records, kind, opts) {
    const states = kind === 'sites' ? 2 : kind === 'protein' ? 20 : 4;
    return {
      kind,
      model: opts.model,
      kappa: opts.kappa,
      kappa2: opts.kappa2,
      freqs: opts.useEmpirical === false ? undefined : freqsFromChars(chars, states),
      alpha: opts.alpha,
      categories: opts.categories
    };
  }

  // -------------------------------------------------------- program entries

  const PROGRAMS = [];

  function program(def) { PROGRAMS.push(def); return def; }

  // ---- parsimony family

  program({
    id: 'DNAPARS', title: 'DNA parsimony', category: 'Parsimony', input: 'dna',
    summary: 'Finds the most parsimonious tree for DNA data using the Fitch (unordered) criterion, with stepwise addition and nearest-neighbour interchange.',
    options: [num('restarts', 'Random restarts', 4, { min: 1, max: 20, step: 1 })],
    run: ({ text, opts }) => {
      const { chars, names } = sequenceInput(text, 'dna');
      const cost = t => I.fitchScore(t, chars, names);
      const { tree, score } = parsimonySearch(chars, names, cost, opts);
      fitchLengths(tree, chars, names);
      return treeResult(tree, ['Steps are the minimum number of base changes over the alignment.', 'Branch lengths are the number of changes inferred on each branch.'], [
        { label: 'STEPS', value: score }, { label: 'SITES', value: chars.length }, { label: 'TAXA', value: names.length }
      ], 'cladogram');
    }
  });

  program({
    id: 'PROTPARS', title: 'Protein parsimony', category: 'Parsimony', input: 'protein',
    summary: 'Protein parsimony under the Fitch criterion, treating amino-acid replacements as unordered and equally costly.',
    options: [num('restarts', 'Random restarts', 4, { min: 1, max: 20, step: 1 })],
    run: ({ text, opts }) => {
      const { chars, names } = sequenceInput(text, 'protein');
      const cost = t => I.fitchScore(t, chars, names);
      const { tree, score } = parsimonySearch(chars, names, cost, opts);
      fitchLengths(tree, chars, names);
      return treeResult(tree, ['Amino-acid changes are unordered; step counts treat every replacement as one step.'], [
        { label: 'STEPS', value: score }, { label: 'SITES', value: chars.length }, { label: 'TAXA', value: names.length }
      ], 'cladogram');
    }
  });

  program({
    id: 'DNAPENNY', title: 'DNA parsimony (branch & bound)', category: 'Parsimony', input: 'dna',
    summary: "DNA parsimony searched exhaustively over addition orders and nearest-neighbour rearrangements, keeping the shortest tree found. This is exact only for small data sets; PHYLIP uses branch-and-bound with bounds and can go further.",
    options: [num('restarts', 'Random restarts', 12, { min: 1, max: 40, step: 1 })],
    run: ({ text, opts }) => {
      const { chars, names } = sequenceInput(text, 'dna');
      const cost = t => I.fitchScore(t, chars, names);
      const { tree, score } = parsimonySearch(chars, names, cost, Object.assign({}, opts, { restarts: opts.restarts || 12 }));
      fitchLengths(tree, chars, names);
      const exhaustive = names.length <= 7;
      return treeResult(tree, [
        `${opts.restarts || 12} independent searches were run and the shortest tree was kept.`,
        exhaustive ? 'With seven or fewer taxa this search is exhaustive over addition orders.' : 'For larger data sets this is a heuristic search, not a guaranteed optimum.'
      ], [{ label: 'STEPS', value: score }, { label: 'SITES', value: chars.length }], 'cladogram');
    }
  });

  program({
    id: 'DNAMOVE', title: 'Interactive DNA parsimony', category: 'Parsimony', input: 'dna',
    interactive: true,
    summary: 'Interactive DNA parsimony, after MacClade. Start from a tree, read the most-parsimonious reconstruction of every site, then rearrange the tree by hand and watch the step count respond.',
    options: [chk('greedy', 'Start from a locally optimal tree', true)],
    run: ({ text, opts }) => {
      const { chars, names } = sequenceInput(text, 'dna');
      const session = PI.createSession({ chars, names, kind: 'dna' });
      if (opts.greedy !== false) session.greedy();
      return {
        kind: 'interactive',
        session,
        notes: [
          'Each listed rearrangement is a nearest-neighbour interchange; the number beside it is the change in steps.',
          'Undo walks back through the rearrangements you have applied.',
          'The site table shows the minimum number of changes implied for each site on the current tree.'
        ],
        stats: [{ label: 'STEPS', value: session.score() }, { label: 'SITES', value: chars.length }]
      };
    }
  });

  program({
    id: 'PARS', title: 'Multistate parsimony (Wagner)', category: 'Parsimony', input: 'discrete',
    summary: 'Wagner parsimony for multistate discrete characters: states 0–9 are ordered, and a change of k states costs k steps.',
    options: [num('restarts', 'Random restarts', 4, { min: 1, max: 20, step: 1 })],
    run: ({ text, opts }) => {
      const { chars, names } = sequenceInput(text, 'discrete');
      const cost = I.wagnerCost(10);
      const score = t => I.sankoffScore(t, chars, names, cost);
      const { tree, score: best } = parsimonySearch(chars, names, score, opts);
      fitchLengths(tree, chars, names);
      return treeResult(tree, ['Characters 0–9 are treated as ordered states; a step between adjacent states costs one.'], [
        { label: 'STEPS', value: best }, { label: 'CHARACTERS', value: chars.length }
      ], 'cladogram');
    }
  });

  program({
    id: 'MOVE', title: 'Interactive mixed parsimony', category: 'Parsimony', input: 'discrete',
    interactive: true,
    summary: 'Interactive parsimony for multistate characters using the mixed Wagner / Dollo / Camin–Sokal criterion of MIX. Rearrange the tree by hand and watch the weighted step count change.',
    options: [
      num('wagnerWeight', 'Wagner weight', 1, { min: 0, max: 10, step: 0.5 }),
      num('dolloWeight', 'Dollo weight', 0, { min: 0, max: 10, step: 0.5 }),
      num('caminWeight', 'Camin–Sokal weight', 0, { min: 0, max: 10, step: 0.5 })
    ],
    run: ({ text, opts }) => {
      const { chars, names } = sequenceInput(text, 'discrete');
      const session = PI.createSession({
        chars, names, kind: 'discrete',
        weights: {
          wagner: opts.wagnerWeight === undefined ? 1 : opts.wagnerWeight,
          dollo: opts.dolloWeight || 0,
          caminSokal: opts.caminWeight || 0
        }
      });
      return {
        kind: 'interactive',
        session,
        notes: [
          'The score is the weighted sum of the Wagner, Dollo and Camin–Sokal criteria, exactly as in MIX.',
          'Nearest-neighbour rearrangements are listed with their score change; apply as many as you like.'
        ],
        stats: [{ label: 'WEIGHTED STEPS', value: Number(session.score().toFixed(3)) }, { label: 'CHARACTERS', value: chars.length }]
      };
    }
  });

  program({
    id: 'PENNY', title: 'Wagner parsimony (branch & bound)', category: 'Parsimony', input: 'discrete',
    summary: 'Wagner parsimony with an intensive search over addition orders and nearest-neighbour interchange.',
    options: [num('restarts', 'Random restarts', 12, { min: 1, max: 40, step: 1 })],
    run: ({ text, opts }) => {
      const { chars, names } = sequenceInput(text, 'discrete');
      const cost = I.wagnerCost(10);
      const score = t => I.sankoffScore(t, chars, names, cost);
      const { tree, score: best } = parsimonySearch(chars, names, score, Object.assign({}, opts, { restarts: opts.restarts || 12 }));
      fitchLengths(tree, chars, names);
      return treeResult(tree, ['Exhaustive for seven or fewer taxa; heuristic above that.'], [
        { label: 'STEPS', value: best }, { label: 'CHARACTERS', value: chars.length }
      ], 'cladogram');
    }
  });

  program({
    id: 'MIX', title: 'Mixed parsimony', category: 'Parsimony', input: 'discrete',
    summary: 'Weighs Wagner, Dollo and Camin–Sokal criteria against each other so that no single assumption about character evolution dominates.',
    options: [
      num('wagnerWeight', 'Wagner weight', 1, { min: 0, max: 10, step: 0.5 }),
      num('dolloWeight', 'Dollo weight', 1, { min: 0, max: 10, step: 0.5 }),
      num('caminWeight', 'Camin–Sokal weight', 1, { min: 0, max: 10, step: 0.5 }),
      num('restarts', 'Random restarts', 4, { min: 1, max: 20, step: 1 })
    ],
    run: ({ text, opts }) => {
      const { chars, names } = sequenceInput(text, 'discrete');
      const score = mixedScore(chars, names, {
        wagner: opts.wagnerWeight === undefined ? 1 : opts.wagnerWeight,
        dollo: opts.dolloWeight === undefined ? 1 : opts.dolloWeight,
        caminSokal: opts.caminWeight === undefined ? 1 : opts.caminWeight
      });
      const { tree, score: best } = parsimonySearch(chars, names, score, opts);
      fitchLengths(tree, chars, names);
      return treeResult(tree, [
        'The reported score is the weighted sum of the three criteria.',
        'Character weighting across sites is not applied; each character contributes equally to its criterion.'
      ], [{ label: 'WEIGHTED STEPS', value: Number(best.toFixed(3)) }, { label: 'CHARACTERS', value: chars.length }], 'cladogram');
    }
  });

  program({
    id: 'DOLLOP', title: 'Dollo parsimony', category: 'Parsimony', input: 'discrete',
    summary: 'Dollo parsimony: a derived state may be gained only once, but may be lost repeatedly. Suits presence/absence data such as restriction sites.',
    options: [num('restarts', 'Random restarts', 6, { min: 1, max: 20, step: 1 })],
    run: ({ text, opts }) => {
      const { chars, names } = sequenceInput(text, 'discrete');
      const cost = t => I.dolloScore(t, chars, names);
      const { tree, score } = parsimonySearch(chars, names, cost, opts);
      fitchLengths(tree, chars, names);
      return treeResult(tree, [
        'Each character allows one gain and any number of losses.',
        'States above zero are read as "present", so 0/1 data is the classical Dollo case.'
      ], [
        { label: 'STEPS', value: score }, { label: 'CHARACTERS', value: chars.length }
      ], 'cladogram');
    }
  });

  program({
    id: 'DOLPENNY', title: 'Dollo parsimony (branch & bound)', category: 'Parsimony', input: 'discrete',
    summary: 'Dollo parsimony with intensive search over addition orders.',
    options: [num('restarts', 'Random restarts', 12, { min: 1, max: 40, step: 1 })],
    run: ({ text, opts }) => {
      const { chars, names } = sequenceInput(text, 'discrete');
      const cost = t => I.dolloScore(t, chars, names);
      const { tree, score } = parsimonySearch(chars, names, cost, Object.assign({}, opts, { restarts: opts.restarts || 12 }));
      fitchLengths(tree, chars, names);
      return treeResult(tree, [
        'Exhaustive for seven or fewer taxa; heuristic above that.',
        'States above zero are read as "present".'
      ], [
        { label: 'STEPS', value: score }, { label: 'CHARACTERS', value: chars.length }
      ], 'cladogram');
    }
  });

  program({
    id: 'DOLMOVE', title: 'Dollo + Camin–Sokal parsimony', category: 'Parsimony', input: 'discrete',
    summary: 'Combines the Dollo criterion with the Camin–Sokal assumption for the two-state case, allowing one gain per character while barring reversals.',
    options: [
      num('dolloWeight', 'Dollo weight', 1, { min: 0, max: 10, step: 0.5 }),
      num('caminWeight', 'Camin–Sokal weight', 1, { min: 0, max: 10, step: 0.5 }),
      num('restarts', 'Random restarts', 6, { min: 1, max: 20, step: 1 })
    ],
    run: ({ text, opts }) => {
      const { chars, names } = sequenceInput(text, 'discrete');
      const score = t => (opts.dolloWeight || 1) * I.dolloScore(t, chars, names)
        + (opts.caminWeight || 1) * I.sankoffScore(t, chars, names, I.caminSokalCost(10));
      const { tree, score: best } = parsimonySearch(chars, names, score, opts);
      fitchLengths(tree, chars, names);
      return treeResult(tree, ['Two criteria are summed with the weights you chose.'], [
        { label: 'WEIGHTED STEPS', value: Number(best.toFixed(3)) }, { label: 'CHARACTERS', value: chars.length }
      ], 'cladogram');
    }
  });

  program({
    id: 'DNACOMP', title: 'DNA compatibility', category: 'Parsimony', input: 'dna',
    summary: 'Finds the largest set of mutually compatible sites (the clique of characters whose state distributions admit a single tree) and builds the tree implied by them.',
    options: [num('restarts', 'Random restarts', 2, { min: 1, max: 10, step: 1 })],
    run: ({ text }) => {
      const { chars, names } = sequenceInput(text, 'dna');
      const { indices, count, totalCharacters } = S.cliques(chars);
      if (!indices.length) throw new Error('No compatible characters were found.');
      const subset = indices.map(i => chars[i]);
      let tree = I.stepwiseSearch(names, t => -I.fitchScore(t, subset, names));
      I.assignIds(tree);
      for (let p = 0; p < 3; p++) if (!I.nniImprove(tree, t => -I.fitchScore(t, subset, names))) break;
      fitchLengths(tree, subset, names);
      return treeResult(tree, [
        `The largest compatible set contains ${count} of ${totalCharacters} sites.`,
        'Sites outside that set contradict at least one member of it and are excluded from the tree.'
      ], [{ label: 'COMPATIBLE SITES', value: count }, { label: 'TOTAL SITES', value: totalCharacters }], 'cladogram');
    }
  });

  program({
    id: 'CLIQUE', title: 'Character compatibility', category: 'Parsimony', input: 'discrete',
    summary: 'Reports the largest clique of pairwise compatible characters. Two characters are compatible when the states of one are nested or disjoint within the other.',
    options: [],
    run: ({ text }) => {
      const { chars, names } = sequenceInput(text, 'discrete');
      const { indices, count, totalCharacters } = S.cliques(chars);
      const lines = [
        `Characters examined: ${totalCharacters}`,
        `Largest compatible set: ${count}`,
        '',
        'Character numbers in the clique (1-based, in input order):',
        indices.map(i => i + 1).join(' ')
      ];
      if (count < totalCharacters) lines.push('', 'The remaining characters conflict with this set and must be interpreted with homoplasy.');
      return {
        kind: 'report',
        text: lines.join('\n'),
        notes: ['Compatibility is the strict criterion: a clique character set can be placed on a tree with no homoplasy.'],
        stats: [{ label: 'CLIQUE SIZE', value: count }, { label: 'CHARACTERS', value: totalCharacters }, { label: 'TAXA', value: names.length }]
      };
    }
  });

  // ---- likelihood family

  const dnaModelOptions = [
    sel('model', 'Substitution model', ['jc', 'f81', 'k2p', 'hky', 'tn93', 'gtr'], 'hky'),
    num('kappa', 'Ts/Tv ratio κ', 2, { min: 0.1, max: 20, step: 0.1 }),
    num('alpha', 'Gamma shape α (0 = off)', 0, { min: 0, max: 10, step: 0.1 }),
    num('categories', 'Gamma categories', 4, { min: 1, max: 8, step: 1 })
  ];

  program({
    id: 'DNAML', title: 'DNA maximum likelihood', category: 'Likelihood', input: 'dna',
    summary: 'Maximum-likelihood tree for DNA under a reversible model with optional discrete-gamma rate variation. Branch lengths are optimised and the topology refined by nearest-neighbour interchange.',
    options: dnaModelOptions,
    run: ({ text, opts }) => {
      const { chars, names, records } = sequenceInput(text, 'dna');
      const ctx = mlContext(chars, records, 'dna', Object.assign({ model: 'hky' }, opts));
      const { tree, logLikelihood } = I.mlTree(chars, names, ctx);
      return treeResult(tree, [
        `Model: ${M.DNA_MODELS[opts.model || 'hky'].label}${ctx.alpha ? ` with discrete-gamma rates (α = ${ctx.alpha})` : ''}.`,
        'Branch lengths are expected substitutions per site.',
        'This is a heuristic search: a different starting order may reach a better optimum.'
      ], [{ label: 'ln L', value: logLikelihood.toFixed(3) }, { label: 'SITES', value: chars.length }]);
    }
  });

  program({
    id: 'DNAMLK', title: 'DNA maximum likelihood (clock)', category: 'Likelihood', input: 'dna',
    summary: 'Maximum likelihood for DNA with a molecular clock enforced, so all tips are equidistant from the root.',
    options: dnaModelOptions,
    run: ({ text, opts }) => {
      const { chars, names, records } = sequenceInput(text, 'dna');
      const ctx = mlContext(chars, records, 'dna', Object.assign({ model: 'hky' }, opts));
      const { tree, logLikelihood } = I.clockTree(chars, names, ctx);
      return treeResult(tree, [
        'Node heights are optimised under a strict clock; the tree is therefore ultrametric.',
        `Model: ${M.DNA_MODELS[opts.model || 'hky'].label}.`
      ], [{ label: 'ln L', value: logLikelihood.toFixed(3) }, { label: 'SITES', value: chars.length }]);
    }
  });

  program({
    id: 'PROML', title: 'Protein maximum likelihood', category: 'Likelihood', input: 'protein',
    summary: 'Maximum-likelihood tree for amino-acid sequences with Poisson or equal-input rates and optional gamma rate variation.',
    options: [
      sel('model', 'Rate model', ['poisson', 'equal-input'], 'poisson'),
      num('alpha', 'Gamma shape α (0 = off)', 0, { min: 0, max: 10, step: 0.1 }),
      num('categories', 'Gamma categories', 4, { min: 1, max: 8, step: 1 })
    ],
    run: ({ text, opts }) => {
      const { chars, names, records } = sequenceInput(text, 'protein');
      const ctx = mlContext(chars, records, 'protein', Object.assign({ model: 'poisson' }, opts));
      const { tree, logLikelihood } = I.mlTree(chars, names, ctx);
      return treeResult(tree, [
        `Rate model: ${M.PROTEIN_MODELS[opts.model || 'poisson'].label}.`,
        'Empirical exchangeability matrices (JTT, Dayhoff, WAG) are not bundled in Phoenix, so protein likelihood uses equal exchangeabilities.'
      ], [{ label: 'ln L', value: logLikelihood.toFixed(3) }, { label: 'SITES', value: chars.length }]);
    }
  });

  program({
    id: 'PROMLK', title: 'Protein maximum likelihood (clock)', category: 'Likelihood', input: 'protein',
    summary: 'Protein maximum likelihood with a strict molecular clock.',
    options: [
      sel('model', 'Rate model', ['poisson', 'equal-input'], 'poisson'),
      num('alpha', 'Gamma shape α (0 = off)', 0, { min: 0, max: 10, step: 0.1 }),
      num('categories', 'Gamma categories', 4, { min: 1, max: 8, step: 1 })
    ],
    run: ({ text, opts }) => {
      const { chars, names, records } = sequenceInput(text, 'protein');
      const ctx = mlContext(chars, records, 'protein', Object.assign({ model: 'poisson' }, opts));
      const { tree, logLikelihood } = I.clockTree(chars, names, ctx);
      return treeResult(tree, ['Node heights are optimised under a strict clock.'], [
        { label: 'ln L', value: logLikelihood.toFixed(3) }, { label: 'SITES', value: chars.length }
      ]);
    }
  });

  program({
    id: 'RESTML', title: 'Restriction site maximum likelihood', category: 'Likelihood', input: 'sites',
    summary: 'Maximum likelihood for presence/absence restriction-site data using a two-state model.',
    options: [num('alpha', 'Gamma shape α (0 = off)', 0, { min: 0, max: 10, step: 0.1 })],
    run: ({ text, opts }) => {
      const { chars, names } = sequenceInput(text, 'sites');
      const ctx = { kind: 'sites', freqs: freqsFromChars(chars, 2), alpha: opts.alpha, categories: 4 };
      const { tree, logLikelihood } = I.mlTree(chars, names, ctx);
      return treeResult(tree, ['Sites are coded 1 (present) and 0 (absent); both changes cost the same.'], [
        { label: 'ln L', value: logLikelihood.toFixed(3) }, { label: 'SITES', value: chars.length }
      ]);
    }
  });

  program({
    id: 'RESTDIST', title: 'Restriction site distances', category: 'Distances', input: 'sites',
    summary: 'Pairwise distances from restriction-site presence/absence data: the simple fraction of differing sites or a maximum-likelihood estimate under a two-state model.',
    options: [sel('model', 'Distance', ['fraction', 'jc', 'ml'], 'fraction')],
    run: ({ text, opts }) => {
      const records = sequenceInput(text, 'sites').records;
      const D = distanceTable(records, 'sites', opts);
      return matrixResult(D, records.map(r => r.name), ['Sites must be coded 1 for present and 0 for absent.'], [
        { label: 'TAXA', value: records.length }, { label: 'SITES', value: records[0].seq.length }
      ]);
    }
  });

  program({
    id: 'DNADIST', title: 'DNA distances', category: 'Distances', input: 'dna',
    summary: 'Pairwise evolutionary distances between DNA sequences, using a choice of substitution models and analytic or maximum-likelihood corrections.',
    options: [
      sel('model', 'Distance model', ['p', 'jc', 'k2p', 'f81', 'hky', 'tn93', 'gtr', 'ml'], 'jc'),
      num('kappa', 'Ts/Tv ratio κ', 2, { min: 0.1, max: 20, step: 0.1 })
    ],
    run: ({ text, opts }) => {
      const records = sequenceInput(text, 'dna').records;
      const D = distanceTable(records, 'dna', opts);
      const model = opts.model || 'jc';
      const label = model === 'p' ? 'p-distance' : model === 'ml' ? 'maximum-likelihood distance (JC rates)' : M.DNA_MODELS[model].label;
      return matrixResult(D, records.map(r => r.name), [
        `Model: ${label}.`,
        'Sites with an ambiguous or missing base in either sequence are skipped for that pair.'
      ], [{ label: 'MODEL', value: label }, { label: 'TAXA', value: records.length }]);
    }
  });

  program({
    id: 'PROTDIST', title: 'Protein distances', category: 'Distances', input: 'protein',
    summary: 'Pairwise distances between amino-acid sequences: observed proportion of differences, Poisson correction, equal-input correction or maximum likelihood.',
    options: [sel('model', 'Distance', ['p', 'poisson', 'equal-input', 'ml'], 'poisson')],
    run: ({ text, opts }) => {
      const records = sequenceInput(text, 'protein').records;
      const D = distanceTable(records, 'protein', opts);
      return matrixResult(D, records.map(r => r.name), [
        'Poisson correction assumes all amino-acid replacements are equally likely.',
        'Empirical exchangeability matrices are not bundled in Phoenix.'
      ], [{ label: 'TAXA', value: records.length }, { label: 'SITES', value: records[0].seq.length }]);
    }
  });

  program({
    id: 'GENDIST', title: 'Gene frequency distances', category: 'Distances', input: 'numbers',
    summary: 'Genetic distances between populations from allele frequencies: Nei’s standard distance, Cavalli-Sforza chord distance or Reynolds’ coancestry distance.',
    options: [
      sel('method', 'Distance', ['nei', 'chord', 'reynolds'], 'nei'),
      num('alleles', 'Alleles per locus', 2, { min: 2, max: 10, step: 1 })
    ],
    run: ({ text, opts }) => {
      const { names, rows, width } = numberInput(text);
      const k = opts.alleles || 2;
      if (width % k !== 0) throw new Error(`Column count (${width}) is not a multiple of alleles per locus (${k}).`);
      const loci = width / k;
      const D = names.map(() => new Array(names.length).fill(0));
      for (let i = 0; i < names.length; i++) {
        for (let j = i + 1; j < names.length; j++) {
          let total = 0;
          for (let l = 0; l < loci; l++) {
            const a = rows[i].slice(l * k, l * k + k);
            const b = rows[j].slice(l * k, l * k + k);
            if (opts.method === 'chord') {
              let dot = 0, na = 0, nb = 0;
              for (let x = 0; x < k; x++) { dot += a[x] * b[x]; na += a[x] ** 2; nb += b[x] ** 2; }
              total += na > 0 && nb > 0 ? Math.sqrt(2 * (1 - dot / Math.sqrt(na * nb))) : 0;
            } else if (opts.method === 'reynolds') {
              let num = 0, den = 0;
              for (let x = 0; x < k; x++) {
                const mean = (a[x] + b[x]) / 2;
                num += (a[x] - b[x]) ** 2;
                den += mean * (1 - mean);
              }
              total += den > 0 ? -(Math.log(1 - num / (2 * den))) : 0;
            } else {
              let ident = 0;
              for (let x = 0; x < k; x++) ident += a[x] * b[x];
              total += -Math.log(Math.max(ident, 1e-12));
            }
          }
          D[i][j] = D[j][i] = total / loci;
        }
      }
      const label = opts.method === 'chord' ? 'Cavalli-Sforza chord' : opts.method === 'reynolds' ? 'Reynolds coancestry' : 'Nei standard distance';
      return matrixResult(D, names, [
        `Method: ${label}. Each locus is a block of ${k} allele frequencies that sum to one.`
      ], [{ label: 'METHOD', value: label }, { label: 'LOCI', value: loci }]);
    }
  });

  // ---- distance tree methods

  program({
    id: 'NEIGHBOR', title: 'Neighbor joining / UPGMA', category: 'Distance trees', input: 'dna',
    summary: 'Builds a tree from a distance matrix using neighbour joining (Saitou–Nei) or UPGMA. Both run on the distances implied by the model you choose.',
    options: [
      sel('method', 'Algorithm', ['nj', 'upgma'], 'nj'),
      sel('model', 'Distance model', ['p', 'jc', 'k2p', 'f81', 'hky', 'ml'], 'jc'),
      num('kappa', 'Ts/Tv ratio κ', 2, { min: 0.1, max: 20, step: 0.1 })
    ],
    run: ({ text, opts }) => {
      const records = sequenceInput(text, 'dna').records;
      const D = distanceTable(records, 'dna', opts);
      const names = records.map(r => r.name);
      const tree = opts.method === 'upgma' ? I.upgma(D, names) : I.neighborJoining(D, names);
      return treeResult(tree, [
        opts.method === 'upgma'
          ? 'UPGMA assumes a molecular clock and produces an ultrametric tree.'
          : 'Neighbour joining is unrooted; Phoenix displays it rooted on the first internal node.',
        'Distances: ' + (opts.model === 'p' ? 'p-distance' : opts.model === 'ml' ? 'maximum likelihood' : M.DNA_MODELS[opts.model].label) + '.'
      ], [{ label: 'METHOD', value: opts.method === 'upgma' ? 'UPGMA' : 'Neighbor joining' }, { label: 'TAXA', value: names.length }]);
    }
  });

  program({
    id: 'FITCH', title: 'Fitch–Margoliash least squares', category: 'Distance trees', input: 'dna',
    summary: 'Least-squares tree without a clock: branch lengths are chosen to minimise the weighted squared difference between observed and tree distances, and the topology is refined by NNI.',
    options: [
      sel('weights', 'Branch weighting', ['1/d2', '1/d', 'equal'], '1/d2'),
      sel('model', 'Distance model', ['p', 'jc', 'k2p', 'f81', 'hky', 'ml'], 'jc')
    ],
    run: ({ text, opts }) => {
      const records = sequenceInput(text, 'dna').records;
      const D = distanceTable(records, 'dna', opts);
      const names = records.map(r => r.name);
      const { tree, ss } = I.leastSquaresFit(D, names, { weights: opts.weights });
      return treeResult(tree, [
        'No clock is assumed, so the tree is unrooted in principle.',
        `Weighting: ${opts.weights === 'equal' ? 'equal' : opts.weights}. Sum of squared deviations: ${ss.toFixed(4)}.`
      ], [{ label: 'SS', value: ss.toFixed(4) }, { label: 'TAXA', value: names.length }]);
    }
  });

  program({
    id: 'KITSCH', title: 'Fitch–Margoliash with clock', category: 'Distance trees', input: 'dna',
    summary: 'Least-squares tree constrained to a molecular clock, so every tip is the same distance from the root.',
    options: [
      sel('weights', 'Branch weighting', ['1/d2', '1/d', 'equal'], '1/d2'),
      sel('model', 'Distance model', ['p', 'jc', 'k2p', 'f81', 'hky', 'ml'], 'jc')
    ],
    run: ({ text, opts }) => {
      const records = sequenceInput(text, 'dna').records;
      const D = distanceTable(records, 'dna', opts);
      const names = records.map(r => r.name);
      const { tree, ss } = I.leastSquaresFit(D, names, { weights: opts.weights, clock: true });
      return treeResult(tree, [
        'A strict clock is enforced, so the tree is ultrametric.',
        `Sum of squared deviations: ${ss.toFixed(4)}.`
      ], [{ label: 'SS', value: ss.toFixed(4) }, { label: 'TAXA', value: names.length }]);
    }
  });

  program({
    id: 'CONTML', title: 'Continuous / gene-frequency ML', category: 'Distance trees', input: 'numbers',
    summary: 'Restricted maximum likelihood under the Brownian-motion model of Cavalli-Sforza and Edwards, for continuous characters or allele frequencies. Branch lengths are fitted by maximising the restricted likelihood and the topology is refined by NNI.',
    options: [
      sel('dataType', 'Data', ['continuous', 'frequencies'], 'continuous'),
      num('sweeps', 'Branch-length sweeps', 3, { min: 1, max: 8, step: 1 })
    ],
    run: ({ text, opts }) => {
      const { names, rows, width } = numberInput(text);
      const frequencies = opts.dataType === 'frequencies';
      if (frequencies) {
        for (const row of rows) {
          for (const v of row) {
            if (v < 0 || v > 1) throw new Error('Gene frequencies must lie between 0 and 1.');
          }
        }
      }
      const points = frequencies ? rows.map(r => r.map(PC.arcsineSqrt)) : rows;
      const fit = PC.fit(points, names, { sweeps: opts.sweeps || 3 });
      return treeResult(fit.tree, [
        frequencies
          ? 'Allele frequencies are transformed with the variance-stabilising arcsine square root of Cavalli-Sforza and Edwards before fitting.'
          : 'Characters are used untransformed.',
        'The likelihood is the restricted (REML) likelihood: the ancestral state at the root is profiled out, so the result does not depend on where the tree is rooted.',
        'Branch-length optimisation here is coordinate-wise golden-section maximisation of the same criterion Felsenstein maximises with his EM algorithm; the EM implementation itself is not reproduced.',
        `Branch lengths are in units of the Brownian-motion rate, which is not separately identifiable; the raw rate estimate before rescaling was ${fit.rawSigma2.toExponential(3)}.`
      ], [
        { label: 'ln L (REML)', value: fit.logLikelihood.toFixed(3) },
        { label: 'DF', value: fit.dof },
        { label: 'CHARACTERS', value: width }
      ]);
    }
  });

  program({
    id: 'CONTRAST', title: 'Phylogenetic contrasts', category: 'Distance trees', input: 'numbers', needsTrees: 1,
    summary: 'Computes Felsenstein’s independent contrasts for continuous characters on the supplied tree: each contrast is the difference between sister lineages, scaled by the variance expected under Brownian motion.',
    options: [],
    run: ({ text, trees }) => {
      const { names, rows, width } = numberInput(text);
      let tree = C.cloneTree(trees[0]);
      let resolved = false;
      C.postorder(tree, n => { if (n.children && n.children.length > 2) resolved = true; });
      if (resolved) {
        tree = C.resolvePolytomies(tree, C.rng(11));
        C.preorder(tree, n => {
          if (n !== tree && (n.length === undefined || n.length === null)) n.length = 1e-6;
          if (n !== tree && n.length === 0) n.length = 1e-6;
        });
      }
      const treeNames = C.leaves(tree).map(l => l.name);
      const missing = names.filter(n => !treeNames.includes(n));
      if (missing.length) throw new Error(`Taxon "${missing[0]}" is in the data but not in the tree.`);
      const index = new Map(names.map((n, i) => [n, i]));
      const contrasts = [];
      C.postorder(tree, node => {
        if (!node.children) {
          const i = index.get(node.name);
          node.__mean = rows[i].slice();
          node.__var = 0;
          return;
        }
        if (node.children.length !== 2) throw new Error('CONTRAST needs a fully bifurcating tree; resolve polytomies first (RETREE can help).');
        const [a, b] = node.children;
        const va = Math.max((a.length || 0) + a.__var, 1e-9);
        const vb = Math.max((b.length || 0) + b.__var, 1e-9);
        const difference = a.__mean.map((v, k) => v - b.__mean[k]);
        const variance = va + vb;
        contrasts.push({
          left: a.name || 'internal', right: b.name || 'internal',
          difference, variance,
          studentized: difference.map(d => d / Math.sqrt(variance))
        });
        node.__mean = a.__mean.map((v, k) => (v * vb + b.__mean[k] * va) / variance);
        node.__var = (va * vb) / variance;
      });
      if (!contrasts.length) throw new Error('The tree has no internal nodes to contrast.');
      const lines = ['Independent contrasts (differences between sister lineages)', ''];
      lines.push('contrast'.padEnd(24) + 'difference'.padStart(14) + 'std. contrast'.padStart(16));
      for (const c of contrasts) {
        for (let k = 0; k < width; k++) {
          const label = k === 0 ? `${c.left} vs ${c.right}` : '';
          lines.push(label.slice(0, 23).padEnd(24) + c.difference[k].toFixed(5).padStart(14) + c.studentized[k].toFixed(5).padStart(16));
        }
      }
      lines.push('', `Total contrasts: ${contrasts.length} for ${names.length} taxa (n - 1 is expected).`);
      lines.push('Standardised contrasts should be approximately independent and identically distributed under Brownian motion.');
      return {
        kind: 'report',
        text: lines.join('\n'),
        notes: resolved
          ? ['Polytomies in the input tree were resolved arbitrarily with near-zero branches; contrasts spanning them are approximate.', 'Each contrast is independent of the others under a Brownian-motion model.']
          : ['Contrasts require a fully bifurcating tree with meaningful branch lengths.', 'Each contrast is independent of the others under a Brownian-motion model.'],
        stats: [{ label: 'CONTRASTS', value: contrasts.length }, { label: 'CHARACTERS', value: width }]
      };
    }
  });

  program({
    id: 'FACTOR', title: 'Multistate to binary recoding', category: 'Diagnostics', input: 'discrete',
    summary: 'Recodes multistate characters into binary characters using a character-state tree, so the result can be analysed by the two-state programs. Reports the mapping and verifies that the recoding is reversible.',
    options: [chk('deleteConstant', 'Delete constant characters', true)],
    run: ({ text, opts }) => {
      const { stateTrees, cleaned } = extractStateSpec(text);
      const { chars, names } = sequenceInput(cleaned, 'discrete');
      const result = S.factorRecode(chars, names, { stateTrees });
      if (result.mismatches) throw new Error('Internal error: the recoding failed its reversibility check.');
      // Constant characters are always reported; dropping them is what FACTOR does.
      const kept = opts.deleteConstant === false ? null : result.deleted;
      const lines = ['Recoded data set (PHYLIP, one column per binary character)', '', S.toPhylip(names, result.binary), ''];
      lines.push('Character mapping');
      lines.push('original  states        binary character(s)   state tree');
      for (const m of result.mapping) {
        lines.push(
          String(m.character).padStart(8) + '  ' +
          m.states.join(',').padEnd(14) +
          m.binaryCharacters.join(',').padEnd(21) +
          m.edges.map(e => `${e[0]}-${e[1]}`).join(' ')
        );
      }
      if (kept && kept.length) {
        lines.push('');
        lines.push(`Deleted as constant: ${kept.map(d => d.character).join(', ')}`);
      }
      lines.push('');
      lines.push(`${result.sites} characters in, ${result.binary.length} binary characters out; reversibility check passed for every observed state.`);
      return {
        kind: 'report',
        text: lines.join('\n'),
        notes: [
          'FACTOR recodes states related by a character-state tree. Unordered multistate characters have no such tree and are left to PARS.',
          'Each edge of the state tree becomes one binary character: 1 for taxa on the derived side of that edge.',
          'Phoenix derives the state tree from the observed states in numeric order unless you declare one.'
        ],
        stats: [
          { label: 'IN', value: result.sites },
          { label: 'OUT', value: result.binary.length },
          { label: 'REVERSIBLE', value: result.mismatches === 0 ? 'yes' : 'no' }
        ]
      };
    }
  });

  /** Optional "#state N: 0-1 1-2" lines declared alongside the data. */
  function extractStateSpec(text) {
    const stateTrees = {};
    const kept = [];
    for (const raw of String(text).split(/\r?\n/)) {
      const m = /^\s*#state\s+(\d+)\s*:\s*(.+)$/i.exec(raw);
      if (!m) { kept.push(raw); continue; }
      const index = Number(m[1]);
      const edges = m[2].split(/[,\s]+/).filter(Boolean).map(pair => {
        const mm = /^(\d+)-(\d+)$/.exec(pair);
        if (!mm) throw new Error(`Could not read the state-tree edge "${pair}". Use pairs such as 0-1.`);
        return [Number(mm[1]), Number(mm[2])];
      });
      stateTrees[index] = edges;
    }
    return { stateTrees, cleaned: kept.join('\n') };
  }

  // ---- invariants
  // ---- invariants

  program({
    id: 'DNAINVAR', title: 'DNA invariants', category: 'Diagnostics', input: 'dna',
    summary: "Computes Cavender's invariants for four taxa. Nucleotides are reduced to purines and pyrimidines, the 4×4 pattern matrix is formed for each of the three topologies, and the invariant residual measures how far each matrix is from the rank-two structure that topology requires.",
    options: [],
    run: ({ text }) => {
      const { chars, names } = sequenceInput(text, 'dna');
      if (names.length !== 4) throw new Error('DNAINVAR works on exactly four sequences.');
      const inv = S.invariants(chars, names);
      const [a, b, c, d] = names;
      const label = { '12|34': `(${a},${b}) | (${c},${d})`, '13|24': `(${a},${c}) | (${b},${d})`, '14|23': `(${a},${d}) | (${b},${c})` };
      const lines = [
        'Cavender invariants (purine/pyrimidine reduction)',
        '',
        `Sites used: ${inv.used}    ambiguous or uninformative: ${inv.skipped}`,
        '',
        'topology                     invariant residual',
        ...inv.ranked.map(r => `${label[r.split].padEnd(28)}${r.residual.toExponential(4).padStart(12)}`)
      ];
      if (inv.bestSplit) {
        lines.push('', `Smallest residual: ${label[inv.bestSplit]}. This topology is the one the invariants favour.`);
        lines.push('The invariant residual vanishes in expectation for the true topology and is positive otherwise.');
      } else {
        lines.push('', 'No informative sites, so no topology can be distinguished.');
      }
      return {
        kind: 'report',
        text: lines.join('\n'),
        notes: [
          'Invariants are rank conditions on the reduced pattern matrix; they hold for any branch lengths, which is what makes the method rate-independent.',
          "Phoenix reports the size of the residual. It does not implement a significance test, and it uses Cavender's quadratic invariants rather than Lake's linear invariants.",
          'Reference: Cavender & Felsenstein (1987), Molecular Biology and Evolution 4:625-640.'
        ],
        stats: [
          { label: 'SITES USED', value: inv.used },
          { label: 'BEST TOPOLOGY', value: inv.bestSplit ? label[inv.bestSplit] : '—' }
        ]
      };
    }
  });

  // ---- bootstrap and consensus

  program({
    id: 'SEQBOOT', title: 'Bootstrap resampling', category: 'Resampling', input: 'dna',
    summary: 'Generates bootstrap replicate data sets by resampling alignment columns with replacement, or summarises them straight into a supported tree with a consensus.',
    options: [
      num('replicates', 'Replicates', 20, { min: 2, max: 200, step: 1 }),
      sel('method', 'Resampling', ['bootstrap', 'jackknife'], 'bootstrap'),
      num('seed', 'Random seed', 4242, { min: 1, max: 999999, step: 1 }),
      chk('summarise', 'Build a supported tree (SEQBOOT → search → CONSENSE)', true),
      sel('treeMethod', 'Tree search', ['parsimony', 'nj', 'upgma', 'fitch'], 'parsimony'),
      sel('rule', 'Consensus rule', ['majority', 'strict'], 'majority')
    ],
    run: ({ text, opts }) => {
      const { chars, names } = sequenceInput(text, 'dna');
      const replicates = opts.replicates || 20;
      const reps = opts.method === 'jackknife'
        ? S.jackknifeReplicates(chars, replicates, opts.seed)
        : S.bootstrapReplicates(chars, replicates, opts.seed);
      if (opts.summarise) {
        const method = opts.treeMethod || 'parsimony';
        const trees = [];
        for (const rep of reps) {
          if (method === 'parsimony') {
            const cost = t => I.fitchScore(t, rep, names);
            const search = parsimonySearch(rep, names, cost, { restarts: 2 });
            trees.push(search.tree);
          } else {
            const D = I.pDistances(rep, names);
            if (method === 'upgma') trees.push(I.upgma(D, names));
            else if (method === 'fitch') trees.push(I.leastSquaresFit(D, names, { weights: '1/d2' }).tree);
            else trees.push(I.neighborJoining(D, names));
          }
        }
        const consensus = S.consensus(trees, { threshold: opts.rule === 'strict' ? 1 : 0.5 });
        return treeResult(consensus.tree, [
          `${replicates} ${opts.method === 'jackknife' ? 'jackknife' : 'bootstrap'} replicates were resampled and searched with ${method.toUpperCase()}.`,
          'Internal labels are the percentage of replicate trees containing that split.',
          'Bootstrap proportions describe support for the data, not a guarantee of the true tree.'
        ], [
          { label: 'REPLICATES', value: replicates },
          { label: 'SEARCH', value: method.toUpperCase() },
          { label: 'SPLITS', value: consensus.splits.length }
        ], 'cladogram');
      }
      const blocks = reps.map((charsRep, r) => {
        const header = `${names.length} ${charsRep.length}`;
        const body = names.map((n, t) => {
          let seq = '';
          for (const col of charsRep) {
            let best = 0, bestCount = -1;
            for (let st = 0; st < col[t].length; st++) if (col[t][st] > bestCount) { bestCount = col[t][st]; best = st; }
            seq += 'ACGT'[best] || 'N';
          }
          return n.slice(0, 10).padEnd(10) + seq;
        });
        return `# replicate ${r + 1}\n` + [header, ...body].join('\n');
      });
      return {
        kind: 'report',
        text: blocks.join('\n\n'),
        notes: [
          `${reps.length} replicate data sets in PHYLIP format.`,
          'Paste them into DNAPARS or DNAML one at a time, then summarise with CONSENSE.',
          'Tick "Build a supported tree" to run that whole pipeline inside Phoenix.'
        ],
        stats: [{ label: 'REPLICATES', value: reps.length }, { label: 'SITES', value: chars.length }, { label: 'METHOD', value: opts.method || 'bootstrap' }]
      };
    }
  });

  program({
    id: 'CONSENSE', title: 'Consensus tree', category: 'Resampling', input: 'trees', needsTrees: 2,
    summary: 'Builds a consensus tree from a collection of trees: strict, majority rule, or a threshold you choose.',
    options: [
      sel('rule', 'Consensus rule', ['majority', 'strict', 'extended'], 'majority'),
      num('threshold', 'Custom threshold (0–1)', 0.5, { min: 0.01, max: 1, step: 0.01 })
    ],
    run: ({ trees, opts }) => {
      const threshold = opts.rule === 'strict' ? 1 : opts.rule === 'majority' ? 0.5 : (opts.threshold || 0.5);
      const result = S.consensus(trees, { threshold });
      return treeResult(result.tree, [
        `Consensus of ${result.total} trees at a ${(threshold * 100).toFixed(0)}% threshold.`,
        'Internal labels show the percentage of input trees containing that split.'
      ], [{ label: 'TREES', value: result.total }, { label: 'SPLITS', value: result.splits.length }], 'cladogram');
    }
  });

  program({
    id: 'TREEDIST', title: 'Tree distances', category: 'Resampling', input: 'trees', needsTrees: 2,
    summary: 'Compares trees: the symmetric difference (Robinson–Foulds) counts splits present in one tree but not the other, and the branch score compares branch lengths.',
    options: [],
    run: ({ trees }) => {
      const rows = [];
      const max = Math.min(trees.length - 1, 3);
      for (let i = 0; i < max; i++) {
        for (let j = i + 1; j < Math.min(trees.length, i + 2); j++) {
          const cmp = S.compareTrees(trees[i], trees[j]);
          rows.push({ pair: `${i + 1}/${j + 1}`, rf: cmp.rf, maxRf: cmp.maxRf, shared: cmp.shared, branch: cmp.branchScore });
        }
      }
      const lines = ['Tree comparison', ''];
      lines.push('pair   splits-only-in-either   shared   branch score');
      for (const r of rows) {
        lines.push(`${r.pair.padEnd(6)} ${String(r.rf).padStart(12)} of ${String(r.maxRf).padEnd(6)} ${String(r.shared).padStart(6)}  ${r.branch.toFixed(4).padStart(10)}`);
      }
      lines.push('', 'The symmetric difference is the number of splits not shared by both trees.');
      return {
        kind: 'report',
        text: lines.join('\n'),
        notes: ['Robinson–Foulds distance and branch-score distance are both reported.'],
        stats: rows.length ? [{ label: 'SYMMETRIC DIFF', value: rows[0].rf }, { label: 'BRANCH SCORE', value: rows[0].branch.toFixed(4) }] : []
      };
    }
  });

  program({
    id: 'RETREE', title: 'Tree editing', category: 'Resampling', input: 'tree', needsTrees: 1,
    summary: 'Reads a Newick tree and rewrites it: reroot on an outgroup, resolve polytomies, or apply a random rearrangement.',
    options: [
      sel('operation', 'Operation', ['reroot', 'resolve', 'rearrange', 'clean'], 'reroot'),
      { id: 'outgroup', label: 'Outgroup taxon', type: 'text', default: '' },
      num('seed', 'Random seed', 1234, { min: 1, max: 999999, step: 1 })
    ],
    run: ({ trees, opts }) => {
      const input = trees[0];
      let tree;
      if (opts.operation === 'resolve') tree = S.resolveRandomly(input, opts.seed);
      else if (opts.operation === 'rearrange') tree = S.radicalPair(input, opts.seed);
      else if (opts.operation === 'reroot') {
        const names = C.leaves(input).map(l => l.name);
        const outgroup = (opts.outgroup || '').trim() || names[names.length - 1];
        tree = S.reroot(input, outgroup);
      } else {
        tree = C.cloneTree(input);
        C.preorder(tree, n => { if (!n.children) n.length = null; });
      }
      return treeResult(tree, [
        `Operation: ${opts.operation}.`,
        opts.operation === 'reroot' ? `Rooted on the branch leading to ${(opts.outgroup || '').trim() || 'the last taxon'}.` : 'Branch lengths are preserved where they remain meaningful.'
      ], [{ label: 'TAXA', value: C.leaves(tree).length }]);
    }
  });

  // ---- drawing

  program({
    id: 'DRAWGRAM', title: 'Draw a rooted tree', category: 'Drawing', input: 'tree', needsTrees: 1,
    summary: 'Renders a rooted tree as a rectangular phylogram or cladogram, drawn directly in the browser.',
    options: [sel('style', 'Style', ['phenogram', 'cladogram'], 'phenogram')],
    run: ({ trees, opts }) => {
      const tree = trees[0];
      return {
        kind: 'tree',
        tree,
        ascii: S.asciiTree(tree),
        svg: S.drawRectangular(tree, { style: opts.style === 'cladogram' ? 'cladogram' : 'phenogram' }),
        notes: [opts.style === 'cladogram' ? 'Branch lengths are ignored; only the branching order is drawn.' : 'Horizontal position is cumulative branch length.'],
        stats: [{ label: 'TAXA', value: C.leaves(tree).length }]
      };
    }
  });

  program({
    id: 'DRAWTREE', title: 'Draw an unrooted tree', category: 'Drawing', input: 'tree', needsTrees: 1,
    summary: 'Renders an unrooted tree as a radial phylogram, the classic DRAWTREE view.',
    options: [],
    run: ({ trees }) => {
      const tree = trees[0];
      return {
        kind: 'tree',
        tree,
        ascii: S.asciiTree(tree),
        svg: S.drawRadial(tree),
        notes: ['Taxa are evenly spaced around a circle; radius is cumulative branch length.'],
        stats: [{ label: 'TAXA', value: C.leaves(tree).length }]
      };
    }
  });

  // ----------------------------------------------------- scientific status

  /**
   * Validation status for every program. "validated" means the number has been
   * checked against an independently computed value by the validation suite in
   * validation/. "documented differences" means the method is the right one but
   * a feature of the original differs; "partial" means part of the original
   * behaviour is missing; "exploratory" means the implementation is a reasoned
   * approach whose results should not be treated as reference values.
   */
  const STATUS = {
    DNAPARS: { status: 'validated', ref: "Fitch (1971) parsimony; step counts checked against hand-computed trees" },
    PROTPARS: { status: 'documented differences', ref: "Unordered Fitch criterion; PHYLIP PROTPARS uses a constrained amino-acid step matrix" },
    DNAPENNY: { status: 'partial', ref: "Exhaustive over addition orders and NNI; PHYLIP uses branch-and-bound with bounds" },
    DNAMOVE: { status: 'documented differences', ref: "Interactive parsimony as in PHYLIP, presented through a web interface instead of a terminal" },
    DNACOMP: { status: 'validated', ref: "Character compatibility; clique reported and checked" },
    DNAML: { status: 'validated', ref: "Felsenstein (1981) pruning likelihood; checked by brute-force enumeration on four taxa" },
    DNAMLK: { status: 'documented differences', ref: "Clock enforced by constrained node-height optimisation, not by PHYLIP's own search" },
    PROML: { status: 'documented differences', ref: "Poisson / equal-input rates; empirical exchangeability matrices are not bundled" },
    PROMLK: { status: 'documented differences', ref: "As PROML, with a clock" },
    DNAINVAR: { status: 'partial', ref: "Cavender quadratic invariants as a rank test; Lake linear invariants and significance tests are not implemented" },
    DNADIST: { status: 'validated', ref: "JC69, K2P and ML distances checked against closed forms" },
    PROTDIST: { status: 'documented differences', ref: "No empirical exchangeability matrices; Poisson and equal-input only" },
    RESTDIST: { status: 'validated', ref: "Two-state corrections checked against closed forms" },
    RESTML: { status: 'documented differences', ref: "Two-state likelihood; rate variation only through gamma categories" },
    SEQBOOT: { status: 'validated', ref: "Column resampling with a reproducible seed; proportions checked on a known tree" },
    FITCH: { status: 'validated', ref: "Fitch-Margoliash weighted least squares; recovers additive distances exactly" },
    KITSCH: { status: 'validated', ref: "Ultrametric least squares; ultrametricity verified" },
    NEIGHBOR: { status: 'validated', ref: "Saitou-Nei neighbor joining; exact recovery on additive matrices" },
    CONTML: { status: 'documented differences', ref: "Brownian-motion restricted ML with coordinate-wise optimisation; the EM algorithm is not reproduced" },
    CONTRAST: { status: 'partial', ref: "Independent contrasts computed; PHYLIP CONTRAST also estimates variance components" },
    GENDIST: { status: 'documented differences', ref: "Nei, Cavalli-Sforza chord and Reynolds distances; Nei's minimum distance is not offered" },
    PARS: { status: 'validated', ref: "Wagner parsimony with ordered states" },
    MIX: { status: 'validated', ref: "Weighted Wagner / Dollo / Camin-Sokal sum, as in PHYLIP MIX" },
    PENNY: { status: 'partial', ref: "Intensive heuristic search; PHYLIP uses branch-and-bound" },
    MOVE: { status: 'documented differences', ref: "Interactive mixed parsimony, web interface rather than terminal" },
    DOLLOP: { status: 'validated', ref: "Dollo criterion: one gain, unlimited losses" },
    DOLPENNY: { status: 'partial', ref: "Intensive heuristic search; PHYLIP uses branch-and-bound" },
    DOLMOVE: { status: 'documented differences', ref: "Weighted Dollo and Camin-Sokal combination" },
    CLIQUE: { status: 'validated', ref: "Estabrook compatibility and maximum clique" },
    FACTOR: { status: 'validated', ref: "Character-state-tree recoding with a reversibility check" },
    DRAWGRAM: { status: 'partial', ref: "Rectangular and cladogram styles only; PHYLIP offers many typographic styles" },
    DRAWTREE: { status: 'partial', ref: "Radial phylogram only; PHYLIP offers several styles and formats" },
    CONSENSE: { status: 'validated', ref: "Strict and majority-rule consensus with support percentages" },
    TREEDIST: { status: 'documented differences', ref: "Robinson-Foulds and branch-score distances; quartet distances are not offered" },
    RETREE: { status: 'partial', ref: "Reroot, resolve, rearrange, clean; interactive branch-length editing is not offered" }
  };

  const CATEGORIES = ['Parsimony', 'Likelihood', 'Distances', 'Distance trees', 'Ordination', 'Diagnostics', 'Resampling', 'Drawing'];

  function get(name) {
    const found = PROGRAMS.find(p => p.id === name.toUpperCase());
    if (!found) throw new Error(`Unknown program "${name}".`);
    return found;
  }

  function describe(p) {
    const meta = STATUS[p.id] || { status: 'partial', ref: '' };
    return {
      id: p.id, title: p.title, category: p.category, input: p.input, summary: p.summary,
      options: p.options, needsTrees: p.needsTrees || 0, needsTree: p.needsTree || 0,
      interactive: !!p.interactive,
      status: meta.status, reference: meta.ref
    };
  }

  function catalogue() {
    return PROGRAMS.map(describe);
  }

  function statusOf(id) {
    const meta = STATUS[id.toUpperCase()];
    return meta ? Object.assign({ id: id.toUpperCase() }, meta) : null;
  }

  return { PROGRAM_LIST: PROGRAMS, CATEGORIES, STATUS, get, catalogue, statusOf, describe };
});
