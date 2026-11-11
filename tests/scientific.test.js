'use strict';
const test = require('node:test');
const assert = require('node:assert');

const C = require('../phoenix-core.js');
const M = require('../phoenix-models.js');
const I = require('../phoenix-infer.js');
const S = require('../phoenix-support.js');
const PC = require('../phoenix-continuous.js');
const PI = require('../phoenix-interactive.js');
const L = require('../phoenix-linalg.js');
const P = require('../phoenix-programs.js');

// ------------------------------------------------------------- linear algebra

test('Cholesky factorisation inverts a positive definite matrix', () => {
  const A = [[4, 2, 1], [2, 5, 3], [1, 3, 6]];
  const Lc = L.cholesky(A);
  assert.ok(Lc);
  const inv = L.inverse(Lc);
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) {
      let sum = 0;
      for (let k = 0; k < 3; k++) sum += A[i][k] * inv[k][j];
      assert.ok(Math.abs(sum - (i === j ? 1 : 0)) < 1e-12);
    }
  }
});

test('Cholesky rejects a non-positive-definite matrix', () => {
  assert.strictEqual(L.cholesky([[1, 2], [2, 1]]), null);
});

// -------------------------------------------------- interactive parsimony

test('interactive session reports consistent scores on apply and undo', () => {
  const records = C.parseFasta('>A\n0011\n>B\n0111\n>C\n1101\n>D\n1111');
  const chars = C.encode(records, C.alphabet('discrete'));
  const names = records.map(r => r.name);
  const session = PI.createSession({ chars, names, kind: 'discrete' });
  const start = session.score();
  const moves = session.moves();
  assert.ok(moves.length > 0, 'rearrangements are enumerated');
  const target = moves[0];
  const result = session.apply(target.id);
  assert.strictEqual(result.ok, true);
  assert.strictEqual(session.score(), target.score);
  session.undo();
  assert.strictEqual(session.score(), start);
});

test('interactive session can always reach a local optimum', () => {
  const records = C.parseFasta('>A\nACGTA\n>B\nACGTC\n>C\nACGGC\n>D\nACGGC');
  const chars = C.encode(records, C.alphabet('dna'));
  const names = records.map(r => r.name);
  const session = PI.createSession({ chars, names, kind: 'dna' });
  const before = session.score();
  session.greedy();
  const after = session.score();
  assert.ok(after <= before, 'greedy search never increases the step count');
  assert.ok(session.moves().every(m => m.delta >= 0), 'no improving rearrangement remains');
});

test('mixed parsimony sessions honour the criterion weights', () => {
  const records = C.parseFasta('>A\n11\n>B\n10\n>C\n01\n>D\n00');
  const chars = C.encode(records, C.alphabet('discrete'));
  const names = records.map(r => r.name);
  const wagnerOnly = PI.createSession({ chars, names, kind: 'discrete', weights: { wagner: 1 } });
  const dolloOnly = PI.createSession({ chars, names, kind: 'discrete', weights: { dollo: 1 } });
  assert.ok(typeof wagnerOnly.score() === 'number');
  assert.ok(dolloOnly.score() >= wagnerOnly.score(), 'Dollo is never cheaper than Wagner on this data');
});

// ------------------------------------------------------------- Brownian ML

test('restricted likelihood is invariant to the root position', () => {
  const points = [[0.1, 0.2], [0.3, 0.35], [0.6, 0.45], [0.75, 0.7]];
  const tree = C.parseNewick('((A:0.3,B:0.3):0.1,(C:0.3,D:0.3):0.1);');
  const base = PC.remlLogLik(PC.treeCovariance(tree).C, points).logLikelihood;
  for (const outgroup of ['A', 'B', 'C', 'D']) {
    const rerooted = PC.remlLogLik(PC.treeCovariance(S.reroot(tree, outgroup)).C, points).logLikelihood;
    assert.ok(Math.abs(base - rerooted) < 1e-8, `rooting on ${outgroup}`);
  }
});

test('restricted likelihood is invariant to a uniform rescaling of branch lengths', () => {
  const points = [[0.1, 0.2], [0.3, 0.35], [0.6, 0.45], [0.75, 0.7]];
  const a = PC.remlLogLik(PC.treeCovariance(C.parseNewick('((A:0.2,B:0.2):0.3,(C:0.2,D:0.2):0.3);')).C, points);
  const b = PC.remlLogLik(PC.treeCovariance(C.parseNewick('((A:2,B:2):3,(C:2,D:2):3);')).C, points);
  assert.ok(Math.abs(a.logLikelihood - b.logLikelihood) < 1e-8);
});

test('CONTML groups taxa that are close in character space', () => {
  const points = [[0.10, 0.20], [0.11, 0.19], [0.40, 0.10], [0.41, 0.11]];
  const fit = PC.fit(points, ['A', 'B', 'C', 'D'], { sweeps: 2 });
  const rf = C.robinsonFoulds(fit.tree, C.parseNewick('((A,B),(C,D));'));
  assert.strictEqual(rf.rf, 0);
});

test('the arcSine transform maps frequencies into a bounded range', () => {
  assert.ok(Math.abs(PC.arcsineSqrt(0.25) - Math.PI / 6) < 1e-12);
  assert.strictEqual(PC.arcsineSqrt(0), 0);
  assert.ok(Math.abs(PC.arcsineSqrt(1) - Math.PI / 2) < 1e-12);
});

// --------------------------------------------------------------- recoding

test('FACTOR recoding is reversible for every state', () => {
  const records = C.parseFasta('>T1\n0012031230\n>T2\n0011031230\n>T3\n0021010230\n>T4\n0021010231');
  const chars = C.encode(records, C.alphabet('discrete'));
  const names = records.map(r => r.name);
  const result = S.factorRecode(chars, names);
  assert.strictEqual(result.mismatches, 0);
  assert.ok(result.binary.length >= 1);
  // Every binary column contains only 0/1/?
  for (const col of result.binary) {
    for (const cell of col) assert.ok('01?'.includes(cell.seq), cell.seq);
  }
});

test('FACTOR deletes constant characters and keeps them out of the output', () => {
  const records = C.parseFasta('>A\n000\n>B\n000\n>C\n011\n>D\n011');
  const chars = C.encode(records, C.alphabet('discrete'));
  const names = records.map(r => r.name);
  const result = S.factorRecode(chars, names);
  assert.strictEqual(result.deleted.length, 1, 'the constant character is dropped');
  assert.strictEqual(result.binary.length, 2, 'each of the two variable characters gives one binary character');
});

test('FACTOR honours an explicit character-state tree', () => {
  const records = C.parseFasta('>A\n000\n>B\n010\n>C\n020\n>D\n020');
  const chars = C.encode(records, C.alphabet('discrete'));
  const names = records.map(r => r.name);
  // Branching state tree: 1 is the ancestor of both 0 and 2.
  const result = S.factorRecode(chars, names, { stateTrees: { 1: [[1, 0], [1, 2]] } });
  assert.strictEqual(result.mapping[0].edges.length, 2);
  assert.strictEqual(result.mismatches, 0);
});

// ------------------------------------------------------------- invariants

test('invariant residuals are smallest for the generating topology', () => {
  // Two identical pairs: the split (A,C)|(B,D) must win.
  const records = C.parseFasta('>A\nACGTACGTACGTACGT\n>C\nACGTACGTACGTACGT\n>B\nACGTACGTACGTACGA\n>D\nACGTACGTACGTACGA');
  const chars = C.encode(records, C.alphabet('dna'));
  const names = records.map(r => r.name);
  const inv = S.invariants(chars, names);
  assert.strictEqual(inv.bestSplit, '12|34');
  assert.ok(inv.results['12|34'].sumSquares < inv.results['13|24'].sumSquares);
});

test('invariants require exactly four taxa', () => {
  const records = C.parseFasta('>A\nACG\n>B\nACG\n>C\nACG');
  const chars = C.encode(records, C.alphabet('dna'));
  assert.throws(() => S.invariants(chars, ['A', 'B', 'C']), /exactly four/i);
});

// --------------------------------------------------------- program mapping

test('DNAMOVE and MOVE are interactive parsimony, not Camin-Sokal', () => {
  const dna = P.get('DNAMOVE');
  const move = P.get('MOVE');
  assert.strictEqual(dna.interactive, true);
  assert.strictEqual(move.interactive, true);
  assert.strictEqual(P.statusOf('DNAMOVE').status, 'documented differences');
  assert.match(P.statusOf('DNAMOVE').ref, /interactive/i);
});

test('FACTOR is recoding and CONTML is Brownian-motion likelihood', () => {
  assert.match(P.get('FACTOR').summary, /binary/i);
  assert.match(P.get('CONTML').summary, /Brownian/i);
  const factor = P.get('FACTOR').run({ text: '>A\n011\n>B\n010\n>C\n001', opts: { deleteConstant: true } });
  assert.strictEqual(factor.kind, 'report');
  assert.match(factor.text, /Recoded data set/);
  const contml = P.get('CONTML').run({ text: '4 2\nA 0.1 0.4\nB 0.12 0.38\nC 0.3 0.2\nD 0.28 0.22', opts: { dataType: 'continuous', sweeps: 2 } });
  assert.strictEqual(contml.kind, 'tree');
  assert.match(contml.notes.join(' '), /restricted/i);
});

test('every program declares a scientific status and a reference', () => {
  for (const program of P.catalogue()) {
    assert.ok(program.status, `${program.id} has a status`);
    assert.ok(['validated', 'documented differences', 'partial', 'exploratory'].includes(program.status), program.status);
    assert.ok(program.reference && program.reference.length > 10, `${program.id} has a reference note`);
  }
});

test('programs marked validated are the ones with numerical checks', () => {
  const validated = P.catalogue().filter(p => p.status === 'validated').map(p => p.id);
  for (const id of ['DNAPARS', 'DNADIST', 'NEIGHBOR', 'FITCH', 'KITSCH', 'CONSENSE', 'FACTOR', 'CLIQUE', 'PARS', 'MIX', 'DOLLOP']) {
    assert.ok(validated.includes(id), `${id} should be validated`);
  }
});
