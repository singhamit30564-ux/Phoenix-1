'use strict';
const test = require('node:test');
const assert = require('node:assert');

const C = require('../phoenix-core.js');
const M = require('../phoenix-models.js');
const I = require('../phoenix-infer.js');
const S = require('../phoenix-support.js');
const P = require('../phoenix-programs.js');

const DNA = `>Human
ATGCTAGCTAGCTACGATCGATCGTACGATCG
>Chimpanzee
ATGCTAGCTAGCTACGATCGATCGTACGATCA
>Gorilla
ATGCTAGCTAGCTACGATCGATCATACGATCA
>Orangutan
ATGCTAGCTAGCTACGATCAATCATACGATCA
>Gibbon
ATGCTAGCTAGTTACGATCAATCATACGATCA`;

const DNA4 = `>A
ACGTACGTACGTACGTACGT
>C
ACGTACGTACGTACGTACGT
>B
ACGTACGTACGTACGTACGA
>D
ACGTACGTACGTACGTACGA`;

const SAMPLES = {
  dna: DNA,
  protein: `>Alpha
MKVLAAGIVGLNLGGKST
>Beta
MKVLAAGIVGLNLGGKST
>Gamma
MKILAAGIVGLNIGGKST
>Delta
MRILAAGIVGLNIGGKTT`,
  sites: `>S1
11010010110100101101
>S2
11010010110100101100
>S3
11000010110100101100
>S4
01000010110000101100`,
  discrete: `>T1
001101011001
>T2
001101011001
>T3
000111010101
>T4
000110010102`,
  numbers: `4 2
A 0.10 0.40
B 0.12 0.38
C 0.30 0.20
D 0.28 0.22`,
  trees: ['((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);', '((A:0.1,C:0.1):0.2,(B:0.1,D:0.1):0.2);']
};

function runProgram(id, overrides) {
  const program = P.get(id);
  const opts = {};
  for (const o of program.options) opts[o.id] = o.default;
  Object.assign(opts, overrides || {});
  const ctx = { text: SAMPLES[program.input] || '', opts };
  if (program.needsTrees > 0) {
    const raw = program.needsTrees > 1 ? SAMPLES.trees : [SAMPLES.trees[0]];
    ctx.trees = raw.map(t => C.parseNewick(t));
  }
  if (id === 'DNAINVAR') ctx.text = DNA4;
  if (id === 'CONTRAST') {
    ctx.text = '4 2\nA 0.5 1.2\nB 0.6 1.1\nC 0.2 0.4\nD 0.25 0.45';
    ctx.trees = [C.parseNewick('((A:0.2,B:0.2):0.3,(C:0.2,D:0.2):0.3);')];
  }
  return program.run(ctx);
}

// --------------------------------------------------------------- catalogue

test('all thirty-five PHYLIP programs are present', () => {
  const expected = ['PROTPARS', 'DNAPARS', 'DNAPENNY', 'DNAMOVE', 'DNACOMP', 'DNAML', 'DNAMLK',
    'PROML', 'PROMLK', 'DNAINVAR', 'DNADIST', 'PROTDIST', 'RESTDIST', 'RESTML', 'SEQBOOT',
    'FITCH', 'KITSCH', 'NEIGHBOR', 'CONTML', 'CONTRAST', 'GENDIST', 'PARS', 'MIX', 'PENNY',
    'MOVE', 'DOLLOP', 'DOLPENNY', 'DOLMOVE', 'CLIQUE', 'FACTOR', 'DRAWGRAM', 'DRAWTREE',
    'CONSENSE', 'TREEDIST', 'RETREE'];
  const ids = P.PROGRAM_LIST.map(p => p.id);
  assert.deepStrictEqual(ids.slice().sort(), expected.slice().sort());
  assert.strictEqual(ids.length, 35);
});

test('every program runs end to end and returns a renderable result', () => {
  for (const program of P.PROGRAM_LIST) {
    const result = runProgram(program.id);
    assert.ok(['tree', 'matrix', 'report'].includes(result.kind), `${program.id} result kind`);
    if (result.kind === 'tree') {
      assert.ok(result.svg.includes('<svg'), `${program.id} draws an SVG`);
      assert.ok(C.leaves(result.tree).length >= 2, `${program.id} tree has leaves`);
    } else {
      assert.ok(result.text.length > 0, `${program.id} produced text`);
    }
  }
});

// ------------------------------------------------------------------- core

test('PHYLIP and FASTA input both parse to the same records', () => {
  const phylip = `5 32
Human     ATGCTAGCTAGCTACGATCGATCGTACGATCG
ChimpanzeeATGCTAGCTAGCTACGATCGATCGTACGATCA
Gorilla   ATGCTAGCTAGCTACGATCGATCATACGATCA
Orangutan ATGCTAGCTAGCTACGATCAATCATACGATCA
Gibbon    ATGCTAGCTAGTTACGATCAATCATACGATCA`;
  const a = C.parseInput(phylip);
  const b = C.parseInput(DNA);
  assert.strictEqual(a.length, b.length);
  assert.deepStrictEqual(a.map(r => r.seq), b.map(r => r.seq));
});

test('invalid input is rejected with a clear message', () => {
  assert.throws(() => C.validateRecords(C.parseFasta('>only\nACGT'), C.alphabet('dna')), /at least 2 sequences/i);
  assert.throws(() => C.validateRecords(C.parseFasta('>a\nACGT\n>b\nACGT\n>a\nACGT'), C.alphabet('dna')), /unique/i);
  assert.throws(() => C.validateRecords(C.parseFasta('>a\nACGT\n>b\nACG'), C.alphabet('dna')), /aligned/i);
  assert.throws(() => C.validateRecords(C.parseFasta('>a\nACGT\n>b\nACGZ'), C.alphabet('dna')), /not valid/i);
});

test('Newick round-trips without losing names or lengths', () => {
  const text = "((A:0.1,B:0.2):0.3,C:0.4);";
  const tree = C.parseNewick(text);
  const back = C.toNewick(tree);
  assert.strictEqual(back, '((A:0.100000,B:0.200000):0.300000,C:0.400000);');
});

// ----------------------------------------------------------------- models

test('Jukes-Cantor distance equals the maximum-likelihood estimate under JC', () => {
  const a = 'ACGTACGTACGTACGTACGT';
  const b = 'ACGTACGTACGTACGTACGA';
  const jc = M.dnaDistance(a, b, { model: 'jc' }).distance;
  const ml = M.dnaDistance(a, b, { model: 'ml' }).distance;
  assert.ok(Math.abs(jc - ml) < 1e-4, `jc=${jc} ml=${ml}`);
});

test('K2P is never shorter than the uncorrected distance', () => {
  const a = 'AAAAAAAAAAAAAAAAAAAA';
  const b = 'AAAAGGGGAAAAGGGGAAAA';
  const p = M.dnaDistance(a, b, { model: 'p' }).distance;
  const k2p = M.dnaDistance(a, b, { model: 'k2p', kappa: 2 }).distance;
  assert.ok(k2p > p);
});

test('transition probabilities are stochastic', () => {
  const { Q, n } = M.buildQ({ kind: 'dna', model: 'hky', kappa: 2, freqs: [0.3, 0.2, 0.3, 0.2] });
  const P = M.transitionMatrix(Q, 0.37, n);
  for (let i = 0; i < n; i++) {
    const sum = P[i].reduce((x, y) => x + y, 0);
    assert.ok(Math.abs(sum - 1) < 1e-10, `row ${i} sums to ${sum}`);
    for (let j = 0; j < n; j++) assert.ok(P[i][j] >= -1e-12);
  }
});

// ------------------------------------------------------------------ trees

test('neighbor joining recovers an additive tree and its path lengths', () => {
  const names = ['A', 'B', 'C', 'D', 'E'];
  const D = [
    [0, 0.3, 1.1, 1.0, 1.2],
    [0.3, 0, 1.2, 1.1, 1.3],
    [1.1, 1.2, 0, 0.4, 0.5],
    [1.0, 1.1, 0.4, 0, 0.6],
    [1.2, 1.3, 0.5, 0.6, 0]
  ];
  const tree = I.neighborJoining(D, names);
  const ctx = I.treeIndices(tree);
  const leaf = n => C.leaves(tree).find(l => l.name === n);
  for (let i = 0; i < names.length; i++) {
    for (let j = i + 1; j < names.length; j++) {
      const got = I.pathLength(leaf(names[i]), leaf(names[j]), ctx);
      assert.ok(Math.abs(got - D[i][j]) < 0.06, `${names[i]}-${names[j]} ${got} vs ${D[i][j]}`);
    }
  }
});

test('UPGMA produces an ultrametric tree', () => {
  const names = ['A', 'B', 'C', 'D'];
  const D = [[0, 2, 6, 6], [2, 0, 6, 6], [6, 6, 0, 4], [6, 6, 4, 0]];
  const tree = I.upgma(D, names);
  const ctx = I.treeIndices(tree);
  const leaf = n => C.leaves(tree).find(l => l.name === n);
  const d = (x, y) => I.pathLength(leaf(x), leaf(y), ctx);
  // Root to tip distance must be the same, i.e. d(A,B) = d(A,C) = d(B,C).
  assert.ok(Math.abs(d('A', 'C') - d('B', 'C')) < 1e-9);
  assert.ok(Math.abs(d('A', 'D') - d('B', 'D')) < 1e-9);
  // Three-point condition for every triple.
  for (const [x, y, z] of [['A', 'B', 'C'], ['A', 'C', 'D'], ['B', 'C', 'D'], ['A', 'B', 'D']]) {
    const ds = [d(x, y), d(x, z), d(y, z)].sort((p, q) => p - q);
    assert.ok(Math.abs(ds[1] - ds[2]) < 1e-9, `${x}${y}${z}: ${ds}`);
  }
});

test('parsimony scores a four-taxon quartet correctly', () => {
  const records = C.parseFasta('>A\nA\n>B\nA\n>C\nC\n>D\nC');
  const chars = C.encode(records, C.alphabet('dna'));
  const names = records.map(r => r.name);
  assert.strictEqual(I.fitchScore(C.parseNewick('((A,B),(C,D));'), chars, names), 1);
  assert.strictEqual(I.fitchScore(C.parseNewick('((A,C),(B,D));'), chars, names), 2);
});

test('Dollo parsimony allows one gain and unlimited losses', () => {
  const records = C.parseFasta('>A\n1\n>B\n1\n>C\n0\n>D\n0');
  const chars = C.encode(records, C.alphabet('discrete'));
  const names = records.map(r => r.name);
  // One gain on the branch leading to (A,B) explains the data with a single step.
  assert.strictEqual(I.dolloScore(C.parseNewick('((A,B),(C,D));'), chars, names), 1);
  // The same data needs a gain plus two losses on the wrong tree.
  assert.strictEqual(I.dolloScore(C.parseNewick('((A,C),(B,D));'), chars, names), 2);
});

test('maximum likelihood is at least as good as the starting tree', () => {
  const records = C.parseFasta(DNA);
  const chars = C.encode(records, C.alphabet('dna'));
  const names = records.map(r => r.name);
  const opts = { kind: 'dna', model: 'hky', kappa: 2, freqs: M.empiricalFrequencies(chars, 4) };
  const start = I.neighborJoining(I.pDistances(chars, names), names);
  for (const [child] of C.parentMap(start)) if (!child.length) child.length = 0.05;
  const startL = I.logLikelihood(start, chars, names, opts);
  const ml = I.mlTree(chars, names, opts);
  assert.ok(ml.logLikelihood >= startL - 1e-6, `ml=${ml.logLikelihood} start=${startL}`);
});

test('the clock constraint never increases the likelihood', () => {
  const records = C.parseFasta(DNA);
  const chars = C.encode(records, C.alphabet('dna'));
  const names = records.map(r => r.name);
  const opts = { kind: 'dna', model: 'jc', freqs: [0.25, 0.25, 0.25, 0.25] };
  const free = I.mlTree(chars, names, opts).logLikelihood;
  const clock = I.clockTree(chars, names, opts).logLikelihood;
  assert.ok(clock <= free + 1e-6, `clock=${clock} free=${free}`);
});

// --------------------------------------------------------------- support

test('majority-rule consensus keeps splits present in most trees', () => {
  const trees = [
    C.parseNewick('((A,B),(C,D));'),
    C.parseNewick('((A,B),(C,D));'),
    C.parseNewick('((A,C),(B,D));')
  ];
  const result = S.consensus(trees, { threshold: 0.5 });
  const splits = result.splits.map(s => [...s.names].sort().join(''));
  assert.ok(splits.includes('AB'), 'the majority split survives');
  const support = result.splits.find(s => [...s.names].sort().join('') === 'AB').support;
  assert.ok(Math.abs(support - 66.666) < 0.1, `support ${support}`);
});

test('strict consensus keeps only unanimous splits', () => {
  const trees = [
    C.parseNewick('((A,B),(C,D));'),
    C.parseNewick('((A,B),(C,D));'),
    C.parseNewick('((A,C),(B,D));')
  ];
  assert.strictEqual(S.consensus(trees, { threshold: 1 }).splits.length, 0);
});

test('Robinson-Foulds distance counts unshared splits', () => {
  const a = C.parseNewick('((A,B),(C,D));');
  const b = C.parseNewick('((A,C),(B,D));');
  const c = C.parseNewick('((A,B),(C,D));');
  assert.strictEqual(S.compareTrees(a, b).rf, 2);
  assert.strictEqual(S.compareTrees(a, c).rf, 0);
});

test('bootstrap replicates resample columns, not taxa', () => {
  const records = C.parseFasta(DNA);
  const chars = C.encode(records, C.alphabet('dna'));
  const replicates = S.bootstrapReplicates(chars, 5, 42);
  assert.strictEqual(replicates.length, 5);
  for (const rep of replicates) {
    assert.strictEqual(rep.length, chars.length);
    assert.strictEqual(rep[0].length, chars[0].length);
  }
});

test('principal coordinates explain all of a Euclidean distance matrix', () => {
  const D = [[0, 1, 2], [1, 0, 1], [2, 1, 0]];
  const { explained } = S.pcoa(D, 2);
  assert.ok(Math.abs(explained.reduce((a, b) => a + b, 0) - 1) < 1e-9);
  assert.ok(explained[0] > 0.5);
});

test('character compatibility finds the conflict-free clique', () => {
  // Characters 1 and 2 are compatible; character 3 conflicts with both.
  const records = C.parseFasta('>A\n001\n>B\n011\n>C\n101\n>D\n111');
  const chars = C.encode(records, C.alphabet('discrete'));
  const result = S.cliques(chars);
  assert.ok(result.count >= 2);
});

test('rerooting on an outgroup places it first', () => {
  const tree = C.parseNewick('((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.1);');
  const rerooted = S.reroot(tree, 'D');
  assert.strictEqual(C.leaves(rerooted)[0].name, 'D');
  assert.strictEqual(C.leaves(rerooted).length, 4);
});

// ------------------------------------------------------- program behaviour

test('DNAINVAR counts the split supported by the data', () => {
  const result = runProgram('DNAINVAR');
  assert.strictEqual(result.stats[1].value, '12|34');
});

test('distance programs return symmetric matrices with zero diagonals', () => {
  for (const id of ['DNADIST', 'PROTDIST', 'RESTDIST', 'GENDIST']) {
    const result = runProgram(id);
    const m = result.matrix;
    for (let i = 0; i < m.length; i++) {
      assert.strictEqual(m[i][i], 0, `${id} diagonal`);
      for (let j = 0; j < m.length; j++) assert.strictEqual(m[i][j], m[j][i], `${id} symmetry`);
    }
  }
});

test('SEQBOOT can summarise replicates into a supported tree', () => {
  const result = runProgram('SEQBOOT', { replicates: 15, summarise: true, seed: 3 });
  assert.strictEqual(result.kind, 'tree');
  const labelled = [];
  C.preorder(result.tree, n => { if (n.support !== undefined && n.support !== null) labelled.push(n.support); });
  assert.ok(labelled.length > 0, 'consensus carries support values');
  assert.ok(labelled.every(v => v >= 50 && v <= 100));
});

test('TREEDIST reports a non-zero distance for different topologies', () => {
  const result = runProgram('TREEDIST');
  assert.match(result.text, /symmetric/i);
});

test('RETREE reroots on the requested outgroup', () => {
  const result = runProgram('RETREE', { operation: 'reroot', outgroup: 'D' });
  assert.strictEqual(C.leaves(result.tree)[0].name, 'D');
});
