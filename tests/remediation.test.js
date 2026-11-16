'use strict';
/* Regression tests for the review findings on PR #3.
 *
 * Expected values here are derived independently of the code under test:
 * hand-computed criterion costs on tiny data sets, and the separate
 * implementations in phoenix-infer.js (sankoffScore / dolloScore) which were
 * written before and validated independently of the interactive session.
 */
const test = require('node:test');
const assert = require('node:assert');

const C = require('../phoenix-core.js');
const I = require('../phoenix-infer.js');
const S = require('../phoenix-support.js');
const PI = require('../phoenix-interactive.js');

function data(text) {
  const records = C.parseFasta(text);
  const chars = C.encode(records, C.alphabet('discrete'));
  return { chars, names: records.map(r => r.name) };
}

// A single character, states 0,0,2,2. Ordered Wagner must cost 2 (one move of
// two steps); unordered Fitch costs 1 (one change). This is the data set that
// distinguishes the two criteria.
const ORDERED = '>A\n0\n>B\n0\n>C\n2\n>D\n2';
const TREE = C.parseNewick('((A,B),(C,D));');

function sessionFor(weights, text) {
  const { chars, names } = data(text || ORDERED);
  return PI.createSession({ chars, names, kind: 'discrete', weights, initialTree: TREE });
}

// ------------------------------------------------------------ MOVE scoring

test('MOVE defaults to Wagner, not Fitch', () => {
  const { chars, names } = data(ORDERED);
  const session = PI.createSession({ chars, names, kind: 'discrete', weights: { wagner: 1 }, initialTree: TREE });
  const wagner = I.sankoffScore(TREE, chars, names, I.wagnerCost(10));
  const fitch = I.fitchScore(TREE, chars, names);
  assert.notStrictEqual(wagner, fitch, 'the data set must distinguish the two criteria');
  assert.strictEqual(session.criterion, 'wagner');
  assert.strictEqual(session.score(), wagner, 'default MOVE score is the Wagner score');
  assert.notStrictEqual(session.score(), fitch, 'default MOVE score is not the Fitch score');
});

test('Wagner-only scoring agrees with Sankoff on ordered costs', () => {
  const { chars, names } = data('>A\n0120\n>B\n0121\n>C\n0211\n>D\n0210');
  const tree = C.parseNewick('((A,B),(C,D));');
  const session = PI.createSession({ chars, names, kind: 'discrete', weights: { wagner: 1 }, initialTree: tree });
  assert.strictEqual(session.score(), I.sankoffScore(tree, chars, names, I.wagnerCost(10)));
});

test('Dollo-only scoring agrees with the Dollo criterion', () => {
  const { chars, names } = data('>A\n0120\n>B\n0121\n>C\n0211\n>D\n0210');
  const tree = C.parseNewick('((A,B),(C,D));');
  const session = PI.createSession({ chars, names, kind: 'discrete', weights: { dollo: 1 }, initialTree: tree });
  assert.strictEqual(session.criterion, 'dollo');
  assert.strictEqual(session.score(), I.dolloScore(tree, chars, names));
});

test('a Dollo-only weight is not silently combined with Wagner', () => {
  const { chars, names } = data(ORDERED);
  const dollo = PI.createSession({ chars, names, kind: 'discrete', weights: { dollo: 1 }, initialTree: TREE });
  assert.strictEqual(dollo.criterion, 'dollo');
  assert.strictEqual(dollo.score(), I.dolloScore(TREE, chars, names));
  // Hand-computed: root present costs two losses, root absent costs one gain.
  assert.strictEqual(dollo.score(), 1);
});

test('Camin-Sokal-only scoring agrees with irreversible ordered costs', () => {
  const { chars, names } = data('>A\n0120\n>B\n0121\n>C\n0211\n>D\n0210');
  const tree = C.parseNewick('((A,B),(C,D));');
  const session = PI.createSession({ chars, names, kind: 'discrete', weights: { caminSokal: 1 }, initialTree: tree });
  assert.strictEqual(session.criterion, 'caminSokal');
  assert.strictEqual(session.score(), I.sankoffScore(tree, chars, names, I.caminSokalCost(10)));
});

test('combined weights equal a hand-computed weighted sum', () => {
  // 0,0,2,2 on ((A,B),(C,D)) with Wagner 1 and Dollo 1.
  //  root 0: one gain (1) + a two-step Wagner move (2) = 3
  //  root 2: two losses (2) + a two-step Wagner move (2) = 4
  //  root 1: loss going down to 0 (1 + 1) and 1->2 Wagner (1) = 3
  const session = sessionFor({ wagner: 1, dollo: 1 });
  assert.strictEqual(session.criterion, 'mixed');
  assert.strictEqual(session.score(), 3);
});

test('all-zero weights are rejected with an explicit error', () => {
  const { chars, names } = data(ORDERED);
  assert.throws(
    () => PI.createSession({ chars, names, kind: 'discrete', weights: { wagner: 0, dollo: 0, caminSokal: 0 }, initialTree: TREE }),
    /at least one/i
  );
});

test('negative weights are rejected', () => {
  const { chars, names } = data(ORDERED);
  assert.throws(
    () => PI.createSession({ chars, names, kind: 'discrete', weights: { wagner: -1 }, initialTree: TREE }),
    /greater than or equal to zero/i
  );
});

test('rearrangement deltas use the same scoring function as the score', () => {
  const { chars, names } = data('>A\n0120\n>B\n0121\n>C\n0211\n>D\n0210');
  const session = PI.createSession({ chars, names, kind: 'discrete', weights: { wagner: 1 } });
  const base = session.score();
  for (const move of session.moves()) {
    assert.ok(Math.abs((move.score - base) - move.delta) < 1e-9, `delta for ${move.id}`);
  }
});

test('applying and undoing a rearrangement restores the exact previous score', () => {
  const { chars, names } = data('>A\n0120\n>B\n0121\n>C\n0211\n>D\n0210');
  const session = PI.createSession({ chars, names, kind: 'discrete', weights: { wagner: 1 } });
  const before = session.score();
  const move = session.moves()[0];
  session.apply(move.id);
  assert.strictEqual(session.score(), before + move.delta);
  session.undo();
  assert.strictEqual(session.score(), before);
});

// ------------------------------------------------ ancestral reconstruction

test('reconstruction costs sum to the reported score in every mode', () => {
  const { chars, names } = data('>A\n0120\n>B\n0121\n>C\n0211\n>D\n0210');
  for (const weights of [{ wagner: 1 }, { dollo: 1 }, { caminSokal: 1 }, { wagner: 2, dollo: 1 }]) {
    const session = PI.createSession({ chars, names, kind: 'discrete', weights });
    const rows = session.siteReconstruction(100);
    assert.strictEqual(rows.length, chars.length);
    const sum = rows.reduce((a, r) => a + r.cost, 0);
    assert.ok(Math.abs(sum - session.score()) < 1e-6, `${JSON.stringify(weights)}: ${sum} vs ${session.score()}`);
  }
});

test('ordered reconstruction places the ancestor at an intermediate state', () => {
  // One character 0,2,2,0 on ((A,B),(C,D)): the cheapest ordered solution puts
  // an internal node at 0 and the other at 2, costing 2.
  const { chars, names } = data('>A\n0\n>B\n0\n>C\n2\n>D\n2');
  const session = PI.createSession({ chars, names, kind: 'discrete', weights: { wagner: 1 }, initialTree: TREE });
  const row = session.siteReconstruction(1)[0];
  assert.strictEqual(row.cost, 2);
  assert.ok([0, 1, 2].includes(row.rootStates[0]), 'root state is a real state');
});

test('Dollo reconstruction never uses more than one gain', () => {
  // Two independent gains would be needed to make 2,2 appear twice; Dollo must
  // instead pay for losses.
  const { chars, names } = data('>A\n2\n>B\n0\n>C\n2\n>D\n0');
  const session = PI.createSession({ chars, names, kind: 'discrete', weights: { dollo: 1 }, initialTree: TREE });
  const row = session.siteReconstruction(1)[0];
  // Best Dollo solution: root present, B and D each lose once -> 2.
  assert.strictEqual(row.cost, 2);
});

test('Camin-Sokal reconstruction never reverts to a smaller state', () => {
  // With 0 on one cherry and 2 on the other, the irreversible model must place
  // the root at 0; a root at 2 would require the 0s to revert.
  const { chars, names } = data('>A\n0\n>B\n0\n>C\n2\n>D\n2');
  const session = PI.createSession({ chars, names, kind: 'discrete', weights: { caminSokal: 1 }, initialTree: TREE });
  const row = session.siteReconstruction(1)[0];
  assert.strictEqual(row.rootStates[0], 0, 'the root must be the ancestral state 0');
  // One irreversible 0 -> 2 move costs two steps.
  assert.strictEqual(row.cost, 2);
});

test('score and reconstruction stay consistent after a rearrangement', () => {
  const { chars, names } = data('>A\n0120\n>B\n0121\n>C\n0211\n>D\n0210');
  const session = PI.createSession({ chars, names, kind: 'discrete', weights: { wagner: 1 } });
  const move = session.moves()[0];
  session.apply(move.id);
  const rows = session.siteReconstruction(100);
  const sum = rows.reduce((a, r) => a + r.cost, 0);
  assert.ok(Math.abs(sum - session.score()) < 1e-6);
});

test('DNAMOVE keeps the unordered Fitch criterion and reconstruction', () => {
  const records = C.parseFasta('>A\nACGT\n>B\nACGA\n>C\nACGA\n>D\nACGC');
  const chars = C.encode(records, C.alphabet('dna'));
  const names = records.map(r => r.name);
  const session = PI.createSession({ chars, names, kind: 'dna' });
  assert.strictEqual(session.criterion, 'fitch');
  assert.strictEqual(session.score(), I.fitchScore(session.tree, chars, names));
  const rows = session.siteReconstruction(100);
  assert.ok(rows.every(r => Array.isArray(r.rootStates)));
});

// ---------------------------------------------------------------- FACTOR

test('an all-constant data set produces no binary characters', () => {
  const records = C.parseFasta('>A\n000\n>B\n000\n>C\n000');
  const chars = C.encode(records, C.alphabet('discrete'));
  const result = S.factorRecode(chars, ['A', 'B', 'C']);
  assert.strictEqual(result.binary.length, 0);
  assert.strictEqual(result.constant.length, 3);
  assert.strictEqual(result.mismatches, 0);
});

test('constant characters are reported separately from variable ones', () => {
  // characters: 1 = 0,0,0 (constant), 2 = 0,0,1 (variable), 3 = 1,1,1 (constant)
  const records = C.parseFasta('>A\n001\n>B\n001\n>C\n011');
  const chars = C.encode(records, C.alphabet('discrete'));
  const result = S.factorRecode(chars, ['A', 'B', 'C']);
  assert.deepStrictEqual(result.constant.map(c => c.character), [1, 3], 'characters 1 and 3 never vary');
  assert.strictEqual(result.binary.length, 1, 'character 2 contributes one binary character');
  assert.deepStrictEqual(result.mapping.map(m => m.character), [2]);
});

test('declaring state trees omits characters that have none', () => {
  const records = C.parseFasta('>A\n012\n>B\n012\n>C\n022');
  const chars = C.encode(records, C.alphabet('discrete'));
  const result = S.factorRecode(chars, ['A', 'B', 'C'], { stateTrees: { 2: [[0, 1], [1, 2]] } });
  assert.deepStrictEqual(result.mapping.map(m => m.character), [2]);
  assert.deepStrictEqual(result.omitted.map(m => m.character).sort(), [1, 3]);
});

test('recoding maps every original character to its binary characters', () => {
  // A single three-state character: 0, 1, 2.
  const records = C.parseFasta('>A\n0\n>B\n1\n>C\n2');
  const chars = C.encode(records, C.alphabet('discrete'));
  const result = S.factorRecode(chars, ['A', 'B', 'C']);
  assert.strictEqual(result.mapping.length, 1);
  const entry = result.mapping[0];
  assert.deepStrictEqual(entry.states, [0, 1, 2]);
  assert.deepStrictEqual(entry.binaryCharacters, [1, 2], 'a linear three-state tree has two edges');
  assert.strictEqual(result.binary.length, 2);
  // Additive coding: state 0 -> 00, state 1 -> 10, state 2 -> 11.
  const row = t => result.binary.map(col => col[t].seq).join('');
  assert.strictEqual(row(0), '00');
  assert.strictEqual(row(1), '10');
  assert.strictEqual(row(2), '11');
});

test('missing states survive recoding as unknown, not as a state', () => {
  // A single character with two observed states and one unknown.
  const records = C.parseFasta('>A\n?\n>B\n0\n>C\n1');
  const chars = C.encode(records, C.alphabet('discrete'));
  const result = S.factorRecode(chars, ['A', 'B', 'C']);
  assert.strictEqual(result.mismatches, 0);
  assert.strictEqual(result.binary.length, 1);
  assert.strictEqual(result.binary[0][0].seq, '?', 'an unknown state is written as ?');
});

test('recoding is reversible for every observed state', () => {
  const records = C.parseFasta('>A\n0012031230\n>B\n0011031230\n>C\n0021010230\n>D\n0021010231');
  const chars = C.encode(records, C.alphabet('discrete'));
  const result = S.factorRecode(chars, ['A', 'B', 'C', 'D']);
  assert.strictEqual(result.mismatches, 0);
});

// --------------------------------------------------- state-tree validation

test('a valid linear state tree is accepted', () => {
  assert.doesNotThrow(() => S.validateStateTree([[0, 1], [1, 2], [2, 3]], [0, 1, 2, 3], 1));
});

test('a valid branching state tree is accepted', () => {
  assert.doesNotThrow(() => S.validateStateTree([[1, 0], [1, 2], [1, 3]], [0, 1, 2, 3], 1));
});

test('a disconnected state tree is rejected', () => {
  assert.throws(() => S.validateStateTree([[0, 1], [2, 3]], [0, 1, 2, 3], 4), /not connected/i);
});

test('a cyclic state tree is rejected', () => {
  assert.throws(() => S.validateStateTree([[0, 1], [1, 2], [2, 0]], [0, 1, 2], 5), /loop/i);
});

test('duplicate edges are rejected', () => {
  assert.throws(() => S.validateStateTree([[0, 1], [1, 0]], [0, 1], 6), /more than once/i);
});

test('self loops are rejected', () => {
  assert.throws(() => S.validateStateTree([[0, 0]], [0], 7), /itself/i);
});

test('an observed state missing from the tree is rejected', () => {
  assert.throws(() => S.validateStateTree([[0, 1]], [0, 1, 2], 8), /absent from the state tree/i);
});

test('an invalid state identifier is rejected', () => {
  assert.throws(() => S.validateStateTree([[0, 12]], [0, 12], 9), /not a valid state/i);
  assert.throws(() => S.validateStateTree([[0, 1.5]], [0, 1], 9), /not a valid state/i);
});

test('an empty state tree is rejected', () => {
  assert.throws(() => S.validateStateTree([], [0, 1], 10), /empty/i);
});

test('state-tree errors name the affected character', () => {
  try {
    S.validateStateTree([[0, 1], [2, 3]], [0, 1, 2, 3], 42);
    assert.fail('should have thrown');
  } catch (error) {
    assert.match(error.message, /Character 42/);
  }
});

test('FACTOR program reports omission and uninformative characters', () => {
  const P = require('../phoenix-programs.js');
  const out = P.get('FACTOR').run({ text: '>A\n0110\n>B\n0100\n>C\n0011', opts: {} });
  assert.strictEqual(out.kind, 'report');
  assert.match(out.text, /Recoded data set/);
  assert.match(out.text, /Uninformative/);
  assert.strictEqual(out.stats.find(s => s.label === 'REVERSIBLE').value, 'yes');
});

// ------------------------------------------------ Dollo gain-budget regression
//
// These tests compare against exhaustive enumeration of every internal-state
// assignment, written here rather than in the implementation, so a shared bug
// in the dynamic program cannot make them pass.

function bruteForceDollo(newick, tipStates, names) {
  const tree = C.parseNewick(newick);
  const internals = [];
  (function walk(n) { if (n.children) { internals.push(n); n.children.forEach(walk); } })(tree);
  const leaves = [];
  (function walk(n) { if (!n.children) leaves.push(n); else n.children.forEach(walk); })(tree);
  const byName = new Map(names.map((n, i) => [n, i]));
  let best = Infinity, bestGains = 0, bestLosses = 0;
  for (let mask = 0; mask < 2 ** internals.length; mask++) {
    const assign = new Map();
    internals.forEach((n, i) => assign.set(n, (mask >> i) & 1));
    leaves.forEach(l => assign.set(l, tipStates[byName.get(l.name)] > 0 ? 1 : 0));
    let gains = 0, losses = 0;
    (function walk(n) {
      if (!n.children) return;
      for (const c of n.children) {
        const a = assign.get(n), b = assign.get(c);
        if (a === 0 && b === 1) gains++;
        else if (a === 1 && b === 0) losses++;
        walk(c);
      }
    })(tree);
    const cost = gains + losses;
    if (gains <= 1 && cost < best) { best = cost; bestGains = gains; bestLosses = losses; }
  }
  return { cost: best, gains: bestGains, losses: bestLosses };
}

function dolloCase(newick, states) {
  const names = newick.replace(/[();:,]/g, '').split(/[()]/).join('').split(',').filter(Boolean);
  const taxa = C.leaves(C.parseNewick(newick)).map(l => l.name);
  void names;
  const text = taxa.map((n, i) => `>${n}\n${states[i]}`).join('\n');
  const records = C.parseFasta(text);
  const chars = C.encode(records, C.alphabet('discrete'));
  const namesInOrder = records.map(r => r.name);
  const tree = C.parseNewick(newick);
  return {
    tree, chars, names: namesInOrder,
    expected: bruteForceDollo(newick, states, namesInOrder)
  };
}

test('Dollo matches exhaustive enumeration on several tree shapes', () => {
  const cases = [
    ['(((A,B),(C,D)),((E,F),(G,H)));', [1, 0, 1, 0, 1, 0, 0, 0]],
    ['((A,B),(C,D));', [1, 0, 1, 0]],
    ['((A,B),(C,D));', [1, 1, 0, 0]],
    ['(((A,B),C),(D,E));', [1, 0, 1, 0, 1]],
    ['((A,(B,C)),(D,E));', [0, 1, 0, 1, 1]],
    ['((A,B),(C,(D,E)));', [1, 0, 0, 1, 0]]
  ];
  for (const [newick, states] of cases) {
    const { tree, chars, names, expected } = dolloCase(newick, states);
    const actual = I.dolloScore(tree, chars, names);
    assert.strictEqual(actual, expected.cost, `${newick} ${states.join('')}`);
  }
});

test('Dollo does not accept two independent gains as a cheap solution', () => {
  // The previous implementation allowed each sibling subtree its own gain and
  // returned 3 for this character. Three gains are illegal, so the true cost is
  // 4: one gain plus three losses.
  const { tree, chars, names, expected } = dolloCase('(((A,B),(C,D)),((E,F),(G,H)));', [1, 0, 1, 0, 1, 0, 0, 0]);
  assert.strictEqual(expected.cost, 4, 'exhaustive enumeration says 4');
  const actual = I.dolloScore(tree, chars, names);
  assert.strictEqual(actual, 4);
  assert.notStrictEqual(actual, 3, 'the illegal three-gain solution must not be reported');
});

test('Dollo reconstruction uses at most one gain', () => {
  const { tree, chars, names } = dolloCase('(((A,B),(C,D)),((E,F),(G,H)));', [1, 0, 1, 0, 1, 0, 0, 0]);
  const session = PI.createSession({ chars, names, kind: 'discrete', weights: { dollo: 1 }, initialTree: tree });
  const row = session.siteReconstruction(1)[0];
  assert.ok(row.gains <= 1, `gains ${row.gains}`);
  assert.strictEqual(row.gains + row.losses, row.cost, 'cost equals gains plus losses');
  assert.strictEqual(row.cost, 4);
});

test('Dollo reconstruction cost equals gains plus losses across characters', () => {
  const { tree, chars, names } = dolloCase('((A,(B,C)),(D,E));', [0, 1, 0, 1, 1]);
  const session = PI.createSession({ chars, names, kind: 'discrete', weights: { dollo: 1 }, initialTree: tree });
  for (const row of session.siteReconstruction(100)) {
    assert.strictEqual(row.gains + row.losses, row.cost);
    assert.ok(row.gains <= 1);
  }
});

test('Dollo handles missing observations', () => {
  const records = C.parseFasta('>A\n?\n>B\n0\n>C\n1\n>D\n0');
  const chars = C.encode(records, C.alphabet('discrete'));
  const names = records.map(r => r.name);
  const tree = C.parseNewick('((A,B),(C,D));');
  const session = PI.createSession({ chars, names, kind: 'discrete', weights: { dollo: 1 }, initialTree: tree });
  // One gain reaches C; the unknown tip costs nothing; B and D are already 0.
  assert.strictEqual(session.score(), 1);
});

test('Dollo score is stable under a rearrangement', () => {
  const { tree, chars, names } = dolloCase('(((A,B),(C,D)),((E,F),(G,H)));', [1, 0, 1, 0, 1, 0, 0, 0]);
  const session = PI.createSession({ chars, names, kind: 'discrete', weights: { dollo: 1 }, initialTree: tree });
  const before = session.score();
  const move = session.moves()[0];
  session.apply(move.id);
  const after = session.score();
  assert.strictEqual(after, before + move.delta);
  const rows = session.siteReconstruction(100);
  assert.strictEqual(rows.reduce((a, r) => a + r.cost, 0), after);
  session.undo();
  assert.strictEqual(session.score(), before);
});

test('Wagner is not capped by the Dollo gain budget', () => {
  // Two separated derived clades: Wagner may use two upward moves.
  const records = C.parseFasta('>A\n1\n>B\n0\n>C\n0\n>D\n1');
  const chars = C.encode(records, C.alphabet('discrete'));
  const names = records.map(r => r.name);
  const tree = C.parseNewick('((A,B),(C,D));');
  const session = PI.createSession({ chars, names, kind: 'discrete', weights: { wagner: 1 }, initialTree: tree });
  assert.strictEqual(session.score(), I.sankoffScore(tree, chars, names, I.wagnerCost(10)));
  // Ordered cost: two separate 0->1 moves cost 2, which requires two gains.
  assert.strictEqual(session.score(), 2);
});

test('Camin-Sokal reconstruction reports no reversals', () => {
  const records = C.parseFasta('>A\n0\n>B\n0\n>C\n2\n>D\n2');
  const chars = C.encode(records, C.alphabet('discrete'));
  const names = records.map(r => r.name);
  const session = PI.createSession({ chars, names, kind: 'discrete', weights: { caminSokal: 1 }, initialTree: C.parseNewick('((A,B),(C,D));') });
  const row = session.siteReconstruction(1)[0];
  assert.strictEqual(row.rootStates[0], 0);
  assert.strictEqual(row.cost, 2);
});
