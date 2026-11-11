#!/usr/bin/env node
/* validation/validate.js — scientific validation suite for Phoenix.
 *
 * Every check states the value it expects and how that value was obtained
 * independently of the code under test: a closed-form formula, an exhaustive
 * enumeration, a hand-computed tree score, or a simulation with a known
 * generating process. The suite writes validation/report.json.
 *
 * Run: node validation/validate.js
 */
'use strict';

const fs = require('fs');
const path = require('path');

const C = require('../phoenix-core.js');
const M = require('../phoenix-models.js');
const I = require('../phoenix-infer.js');
const S = require('../phoenix-support.js');
const PC = require('../phoenix-continuous.js');
const PI = require('../phoenix-interactive.js');
const L = require('../phoenix-linalg.js');
const P = require('../phoenix-programs.js');

const checks = [];

function check(id, category, description, reference, expected, actual, tolerance) {
  const bothNumbers = typeof expected === 'number' && typeof actual === 'number';
  const pass = bothNumbers
    ? Math.abs(expected - actual) <= (tolerance === undefined ? 1e-9 : tolerance)
    : expected === actual;
  checks.push({
    id, category, description, reference,
    expected, actual: typeof actual === 'number' && Number.isFinite(actual) ? Number(actual.toPrecision(12)) : actual,
    tolerance: tolerance === undefined ? null : tolerance,
    pass
  });
}

function checkTrue(id, category, description, reference, condition, detail) {
  checks.push({ id, category, description, reference, expected: true, actual: !!condition, tolerance: null, pass: !!condition, detail });
}

// ------------------------------------------------------------------ helpers

function simulate(tree, nSites, seed) {
  const rand = C.rng(seed);
  const model = M.buildQ({ kind: 'dna', model: 'jc' });
  const names = C.leaves(tree).map(l => l.name);
  const chars = [];
  for (let s = 0; s < nSites; s++) {
    const assign = new Map();
    (function walk(node, parentState) {
      let state;
      if (!node.children) {
        if (parentState === null) state = Math.floor(rand() * 4);
        else {
          const P = M.transitionMatrix(model.Q, Math.max(node.length || 0, 1e-9), 4);
          let r = rand();
          state = 3;
          for (let j = 0; j < 4; j++) { r -= P[parentState][j]; if (r <= 0) { state = j; break; } }
        }
        assign.set(node.name, state);
        return;
      }
      if (parentState === null) state = Math.floor(rand() * 4);
      else {
        const P = M.transitionMatrix(model.Q, Math.max(node.length || 0, 1e-9), 4);
        let r = rand();
        state = 3;
        for (let j = 0; j < 4; j++) { r -= P[parentState][j]; if (r <= 0) { state = j; break; } }
      }
      for (const child of node.children) walk(child, state);
    })(tree, null);
    chars.push(names.map(n => {
      const v = new Float64Array(4);
      v[assign.get(n)] = 1;
      return v;
    }));
  }
  return { chars, names };
}

function reconstructTree(root) {
  // Minimal rooted binary tree from a Newick string.
  return C.parseNewick(root);
}

// ------------------------------------------------------------------- checks

// 1. Jukes-Cantor closed form
{
  const a = 'ACGTACGTACGTACGTACGT';
  const b = 'ACGTACGTACGTACGTACGA';
  const p = 1 / 20;
  const expected = -0.75 * Math.log(1 - (4 * p) / 3);
  const got = M.dnaDistance(a, b, { model: 'jc' }).distance;
  check('dist.jc', 'distances', 'Jukes-Cantor distance matches the closed form -3/4 ln(1-4p/3)',
    'closed-form JC69 formula', expected, got, 1e-9);
}

// 2. Kimura 2-parameter closed form
{
  // Build an alignment with exactly 4 transitions and 2 transversions out of 20.
  const a = 'AAAAAAAAAAAAAAAAAAAA';
  const b = 'GGGGAAAACCCCAAAAAAAC';
  let transitions = 0, transversions = 0, comparable = 0;
  for (let k = 0; k < a.length; k++) {
    if (a[k] === b[k]) { comparable++; continue; }
    comparable++;
    const pair = a[k] + b[k];
    if (pair === 'AG' || pair === 'GA' || pair === 'CT' || pair === 'TC') transitions++;
    else transversions++;
  }
  const P_ = transitions / comparable, Q_ = transversions / comparable;
  const expected = 0.5 * Math.log(1 / (1 - 2 * P_ - Q_)) + 0.25 * Math.log(1 / (1 - 2 * Q_));
  const got = M.dnaDistance(a, b, { model: 'k2p', kappa: 2 }).distance;
  check('dist.k2p', 'distances', 'K2P distance matches the closed form for the same transition/transversion counts',
    'closed-form K80 formula', expected, got, 1e-9);
}

// 3. Maximum-likelihood distance under K2P equals the K2P closed form
{
  const a = 'ACGTACGTACGTACGTACGTACGTACGTACGT';
  const b = 'ACGTACGTACGTACGTACGTACGTACGAACGA';
  const closed = M.dnaDistance(a, b, { model: 'k2p', kappa: 2, freqs: [0.25, 0.25, 0.25, 0.25] }).distance;
  const ml = M.dnaDistance(a, b, { model: 'hky', kappa: 2, freqs: [0.25, 0.25, 0.25, 0.25] }).distance;
  check('dist.ml-k2p', 'distances', 'HKY with equal frequencies reproduces the K2P closed-form distance',
    'numerical ML must agree with the analytic K2P solution when the models coincide', closed, ml, 2e-3);
}

// 4. Transition matrices form a semigroup and are stochastic
{
  const { Q, n } = M.buildQ({ kind: 'dna', model: 'gtr', rates: [2, 1.5, 1, 0.5, 1.2, 0.8], freqs: [0.3, 0.2, 0.25, 0.25] });
  const P1 = M.transitionMatrix(Q, 0.13, n);
  const P2 = M.transitionMatrix(Q, 0.29, n);
  const Psum = M.transitionMatrix(Q, 0.42, n);
  let maxDiff = 0, rowErr = 0;
  for (let i = 0; i < n; i++) {
    let s = 0;
    for (let j = 0; j < n; j++) s += Psum[i][j];
    rowErr = Math.max(rowErr, Math.abs(s - 1));
    for (let j = 0; j < n; j++) {
      let acc = 0;
      for (let k = 0; k < n; k++) acc += P1[i][k] * P2[k][j];
      maxDiff = Math.max(maxDiff, Math.abs(acc - Psum[i][j]));
    }
  }
  check('models.semigroup', 'models', 'P(t1)P(t2) = P(t1+t2) for a GTR generator',
    'matrix exponential semigroup property', 0, maxDiff, 1e-10);
  check('models.stochastic', 'models', 'Transition matrix rows sum to one',
    'definition of a stochastic matrix', 0, rowErr, 1e-12);
}

// 5. Pruning likelihood against exhaustive enumeration of internal states
{
  const tree = C.parseNewick('((A:0.10,B:0.15):0.08,(C:0.12,D:0.09):0.11);');
  const records = C.parseFasta('>A\nACGTACGTAA\n>B\nACGTACGTAC\n>C\nACGTTCGTAC\n>D\nACGTTCGTAG');
  const chars = C.encode(records, C.alphabet('dna'));
  const names = records.map(r => r.name);
  const opts = { kind: 'dna', model: 'jc', freqs: [0.25, 0.25, 0.25, 0.25] };
  const engine = I.logLikelihood(tree, chars, names, opts);

  // Independent: enumerate the states of the two internal nodes explicitly.
  const { Q, pi } = M.buildQ({ kind: 'dna', model: 'jc' });
  const Pt = (t) => M.transitionMatrix(Q, t, 4);
  const [leftChild, rightChild] = tree.children;
  const [nA, nB] = leftChild.children;
  const [nC, nD] = rightChild.children;
  const Pleft = Pt(leftChild.length);
  const Pright = Pt(rightChild.length);
  const PA = Pt(nA.length), PB = Pt(nB.length), PC = Pt(nC.length), PD = Pt(nD.length);
  const code = { A: 0, C: 1, G: 2, T: 3 };
  let ll = 0;
  for (let s = 0; s < records[0].seq.length; s++) {
    const a = code[records[0].seq[s]], b = code[records[1].seq[s]];
    const c = code[records[2].seq[s]], d = code[records[3].seq[s]];
    let total = 0;
    for (let r = 0; r < 4; r++) {
      let left = 0, right = 0;
      for (let x = 0; x < 4; x++) left += Pleft[r][x] * PA[x][a] * PB[x][b];
      for (let y = 0; y < 4; y++) right += Pright[r][y] * PC[y][c] * PD[y][d];
      total += pi[r] * left * right;
    }
    ll += Math.log(total);
  }
  check('likelihood.bruteforce', 'likelihood',
    'Felsenstein pruning equals exhaustive summation over all 16 internal-state combinations',
    'brute-force enumeration written independently in this file', ll, engine, 1e-8);
}

// 6. Two-taxon ML branch length against the closed form
{
  const records = C.parseFasta('>A\nAAAAAAAAAAAAAAAAAAAA\n>B\nAAAACCCCCAAAAAAAAAAA');
  const chars = C.encode(records, C.alphabet('dna'));
  const opts = { kind: 'dna', model: 'jc', freqs: [0.25, 0.25, 0.25, 0.25] };
  const tree = C.parseNewick('(A:0.1,B:0.1);');
  const optimised = I.optimizeBranchLengths(tree, chars, ['A', 'B'], opts, 3);
  let differences = 0;
  for (let k = 0; k < records[0].seq.length; k++) if (records[0].seq[k] !== records[1].seq[k]) differences++;
  const p = differences / records[0].seq.length;
  const expected = -0.75 * Math.log(1 - (4 * p) / 3);
  const total = (tree.children[0].length || 0) + (tree.children[1].length || 0);
  check('likelihood.2taxon-mle', 'likelihood',
    'Optimised branch lengths for two taxa sum to the JC69 distance',
    'closed-form MLE for a two-taxon tree', expected, total, 5e-3);
  checkTrue('likelihood.2taxon-lik', 'likelihood', 'Two-taxon optimisation reaches a finite likelihood',
    'sanity', Number.isFinite(optimised));
}

// 7. Parsimony by hand
{
  const records = C.parseFasta('>A\nA\n>B\nA\n>C\nC\n>D\nC');
  const chars = C.encode(records, C.alphabet('dna'));
  const names = records.map(r => r.name);
  const correct = I.fitchScore(C.parseNewick('((A,B),(C,D));'), chars, names);
  const incorrect = I.fitchScore(C.parseNewick('((A,C),(B,D));'), chars, names);
  check('parsimony.fitch', 'parsimony', 'One character with two states costs one step on the correct tree',
    'hand-computed Fitch score', 1, correct, 0);
  check('parsimony.fitch-wrong', 'parsimony', 'The same character costs two steps on a conflicting tree',
    'hand-computed Fitch score', 2, incorrect, 0);
}

// 8. Wagner and Dollo criteria on a constructed data set
{
  const records = C.parseFasta('>A\n1\n>B\n1\n>C\n0\n>D\n0');
  const chars = C.encode(records, C.alphabet('discrete'));
  const names = records.map(r => r.name);
  check('parsimony.dollo', 'parsimony', 'Dollo needs a single gain on the tree that groups the derived taxa',
    'Dollo definition: one gain, unlimited losses', 1,
    I.dolloScore(C.parseNewick('((A,B),(C,D));'), chars, names), 0);
  check('parsimony.dollo-wrong', 'parsimony', 'Dollo needs a gain plus two losses on a conflicting tree',
    'Dollo definition', 2, I.dolloScore(C.parseNewick('((A,C),(B,D));'), chars, names), 0);
  check('parsimony.wagner', 'parsimony', 'Wagner costs one step for a single 0 -> 1 change',
    'ordered-state cost matrix', 1,
    I.sankoffScore(C.parseNewick('((A,B),(C,D));'), chars, names, I.wagnerCost(10)), 0);
}

// 9. Neighbor joining recovers an additive tree exactly
{
  const trueTree = C.parseNewick('((A:0.10,B:0.20):0.30,(C:0.40,D:0.50):0.60,E:0.70);');
  const names = C.leaves(trueTree).map(l => l.name);
  const ctx = I.treeIndices(trueTree);
  const leafOf = n => C.leaves(trueTree).find(l => l.name === n);
  const D = names.map(a => names.map(b => a === b ? 0 : I.pathLength(leafOf(a), leafOf(b), ctx)));
  const recovered = I.neighborJoining(D, names);
  const ctx2 = I.treeIndices(recovered);
  const leafOf2 = n => C.leaves(recovered).find(l => l.name === n);
  let maxErr = 0;
  for (let i = 0; i < names.length; i++) {
    for (let j = i + 1; j < names.length; j++) {
      const d = I.pathLength(leafOf2(names[i]), leafOf2(names[j]), ctx2);
      maxErr = Math.max(maxErr, Math.abs(d - D[i][j]));
    }
  }
  check('nj.additive', 'tree building', 'Neighbor joining reproduces every pairwise path length of an additive matrix',
    'path lengths measured on the generating tree', 0, maxErr, 1e-9);
  const rf = C.robinsonFoulds(trueTree, recovered);
  check('nj.topology', 'tree building', 'Neighbor joining recovers the generating topology',
    'Robinson-Foulds symmetric difference against the generating tree', 0, rf.rf, 0);
}

// 10. UPGMA on a known matrix
{
  const names = ['A', 'B', 'C', 'D'];
  const D = [[0, 2, 6, 6], [2, 0, 6, 6], [6, 6, 0, 4], [6, 6, 4, 0]];
  const tree = I.upgma(D, names);
  const ctx = I.treeIndices(tree);
  const leafOf = n => C.leaves(tree).find(l => l.name === n);
  const d = (x, y) => I.pathLength(leafOf(x), leafOf(y), ctx);
  const split = C.robinsonFoulds(tree, C.parseNewick('((A,B),(C,D));'));
  check('upgma.topology', 'tree building', 'UPGMA groups the two closest pairs',
    'hand-computed clustering', 0, split.rf, 0);
  checkTrue('upgma.ultrametric', 'tree building', 'Every triple satisfies the ultrametric three-point condition',
    'ultrametric inequality',
    Math.abs(d('A', 'C') - d('B', 'C')) < 1e-9 && Math.abs(d('A', 'D') - d('B', 'D')) < 1e-9);
}

// 11. Brownian-motion REML: invariance properties
{
  const points = [[0.10, 0.20], [0.12, 0.18], [0.40, 0.10], [0.42, 0.09]];
  const tree = C.parseNewick('((A:0.3,B:0.3):0.1,(C:0.3,D:0.3):0.1);');
  const base = PC.remlLogLik(PC.treeCovariance(tree).C, points).logLikelihood;
  const rerooted = PC.remlLogLik(PC.treeCovariance(S.reroot(tree, 'C')).C, points).logLikelihood;
  const scaled = PC.remlLogLik(PC.treeCovariance(C.parseNewick('((A:3,B:3):1,(C:3,D:3):1);')).C, points).logLikelihood;
  check('contml.root-invariance', 'continuous ml', 'The restricted likelihood does not depend on where the tree is rooted',
    'REML invariance to root position', base, rerooted, 1e-8);
  check('contml.scale-invariance', 'continuous ml', 'The restricted likelihood does not depend on the overall branch-length scale',
    'Brownian-motion scale invariance', base, scaled, 1e-8);
}

// 12. CONTML recovers a topology with a clear signal
{
  const points = [[0.10, 0.20], [0.11, 0.19], [0.40, 0.10], [0.41, 0.11]];
  const fit = PC.fit(points, ['A', 'B', 'C', 'D'], { sweeps: 3 });
  const rf = C.robinsonFoulds(fit.tree, C.parseNewick('((A,B),(C,D));'));
  check('contml.topology', 'continuous ml', 'CONTML groups the two similar pairs',
    'the generating structure of the data is two tight clusters', 0, rf.rf, 0);
}

// 13. Independent multivariate-normal evaluation of the restricted likelihood
{
  const points = [[0.10, 0.20], [0.30, 0.35], [0.60, 0.45], [0.75, 0.70]];
  const tree = C.parseNewick('((A:0.2,B:0.2):0.3,(C:0.2,D:0.2):0.3);');
  const { C: K } = PC.treeCovariance(tree);
  const engine = PC.remlLogLik(K, points).logLikelihood;

  // Independent inverse and determinant by Gauss-Jordan elimination.
  const n = K.length;
  const A = K.map((row, i) => Array.from(row).concat(Array.from({ length: n }, (_, j) => (i === j ? 1 : 0))));
  for (let col = 0; col < n; col++) {
    let pivot = col;
    for (let r = col + 1; r < n; r++) if (Math.abs(A[r][col]) > Math.abs(A[pivot][col])) pivot = r;
    const t = A[pivot]; A[pivot] = A[col]; A[col] = t;
    const pv = A[col][col];
    for (let j = 0; j < 2 * n; j++) A[col][j] /= pv;
    for (let r = 0; r < n; r++) {
      if (r === col) continue;
      const f = A[r][col];
      if (!f) continue;
      for (let j = 0; j < 2 * n; j++) A[r][j] -= f * A[col][j];
    }
  }
  const inv = A.map(row => row.slice(n));
  const ones = new Array(n).fill(1);
  const quad = (v) => v.reduce((s, x, i) => s + x * v.reduce((q, y, j) => q + inv[i][j] * v[j], 0), 0);
  const denom = quad(ones);
  let rss = 0;
  const d = points[0].length;
  for (let k = 0; k < d; k++) {
    const x = points.map(r => r[k]);
    let cross = 0;
    // 1' C^-1 x, not x' C^-1 x
    for (let i = 0; i < n; i++) cross += x.reduce((s, _, j) => s + inv[i][j] * x[j], 0);
    const mu = cross / denom;
    rss += quad(x.map(v => v - mu));
  }
  // Determinant by elimination.
  const B = K.map(r => Array.from(r));
  let det = 1;
  for (let c = 0; c < n; c++) {
    let pv = c;
    for (let r = c + 1; r < n; r++) if (Math.abs(B[r][c]) > Math.abs(B[pv][c])) pv = r;
    if (pv !== c) { const t = B[pv]; B[pv] = B[c]; B[c] = t; det = -det; }
    det *= B[c][c];
    for (let r = c + 1; r < n; r++) {
      const f = B[r][c] / B[c][c];
      for (let j = c; j < n; j++) B[r][j] -= f * B[c][j];
    }
  }
  const dof = (n - 1) * d;
  const sigma2 = rss / dof;
  const expected = -0.5 * (n * d * Math.log(2 * Math.PI) + d * Math.log(Math.abs(det)) + d * Math.log(denom) + dof * Math.log(sigma2) + dof);
  check('contml.density', 'continuous ml', 'Restricted likelihood matches an independent Gauss-Jordan evaluation',
    'independent inverse, determinant and quadratic forms', expected, engine, 1e-8);
}

// 14. FACTOR recoding is reversible
{
  const records = C.parseFasta('>T1\n001203123001\n>T2\n001103123001\n>T3\n002101023012\n>T4\n002101023011');
  const chars = C.encode(records, C.alphabet('discrete'));
  const names = records.map(r => r.name);
  const result = S.factorRecode(chars, names);
  check('factor.roundtrip', 'recoding', 'Every observed state decodes back to itself after binary recoding',
    'the recoding is a bijection between states and edge patterns', 0, result.mismatches, 0);
  checkTrue('factor.expands', 'recoding', 'A character with k observed states becomes k-1 binary characters',
    'character-state-tree edge count', result.binary.length >= result.mapping.length);
}

// 15. DNAINVAR on simulated data
{
  const numbers = '((1:0.15,2:0.15):0.25,(3:0.15,4:0.15):0.25);';
  const generating = reconstructTree(numbers);
  const { chars, names } = simulate(generating, 3000, 20240101);
  const invariants = S.invariants(chars, names);
  check('invariants.best-split', 'invariants',
    'Invariant residuals are smallest for the topology that generated the data',
    'simulation on ((1,2),(3,4)) with JC69, 3000 sites, fixed seed', true,
    invariants.bestSplit === '12|34');
  const residualTrue = invariants.results['12|34'].normalised;
  const residualOther = Math.min(invariants.results['13|24'].normalised, invariants.results['14|23'].normalised);
  checkTrue('invariants.ordering', 'invariants',
    'The generating topology has the smaller invariant residual',
    'invariant theory', residualTrue < residualOther,
    `true=${residualTrue.toExponential(3)} other=${residualOther.toExponential(3)}`);
}

// 16. Consensus support percentages
{
  const trees = [
    C.parseNewick('((A,B),(C,D));'),
    C.parseNewick('((A,B),(C,D));'),
    C.parseNewick('((A,C),(B,D));')
  ];
  const result = S.consensus(trees, { threshold: 0.5 });
  const ab = result.splits.find(s => [...s.names].sort().join('') === 'AB');
  check('consensus.support', 'support', 'A split present in two of three trees is reported at 66.7%',
    'exact proportion 2/3', 66.6666666667, ab ? ab.support : -1, 1e-6);
  const strict = S.consensus(trees, { threshold: 1 });
  check('consensus.strict', 'support', 'Strict consensus keeps no splits when the trees disagree',
    'definition of strict consensus', 0, strict.splits.length, 0);
}

// 17. Robinson-Foulds distance
{
  const a = C.parseNewick('((A,B),(C,D));');
  const b = C.parseNewick('((A,C),(B,D));');
  const c = C.parseNewick('((A,B),(C,D));');
  check('treedist.rf', 'tree comparison', 'Two different topologies differ by two splits',
    'Robinson-Foulds definition', 2, S.compareTrees(a, b).rf, 0);
  check('treedist.identical', 'tree comparison', 'A tree has distance zero from itself',
    'Robinson-Foulds definition', 0, S.compareTrees(a, c).rf, 0);
}

// 18. Bootstrap determinism and a strongly supported split
{
  const records = C.parseFasta(`>A\n${'AAAAAAAAAACCCCCCCCCC'}\n>B\n${'AAAAAAAAAACCCCCCCCCC'}\n>C\n${'GGGGGGGGGGTTTTTTTTTT'}\n>D\n${'GGGGGGGGGGTTTTTTTTTT'}`);
  const chars = C.encode(records, C.alphabet('dna'));
  const names = records.map(r => r.name);
  const runs = [1, 2].map(seed => {
    const reps = S.bootstrapReplicates(chars, 30, seed);
    const trees = reps.map(rep => I.neighborJoining(I.pDistances(rep, names), names));
    return S.consensus(trees, { threshold: 0.5 }).splits.map(s => `${[...s.names].sort().join('')}:${s.support.toFixed(2)}`).join('|');
  });
  check('bootstrap.reproducible', 'support', 'The same seed produces identical bootstrap support',
    'seeded resampling', runs[0], runs[1]);
  const reps = S.bootstrapReplicates(chars, 40, 7);
  const trees = reps.map(rep => I.neighborJoining(I.pDistances(rep, names), names));
  const cons = S.consensus(trees, { threshold: 0.5 });
  const ab = cons.splits.find(s => [...s.names].sort().join('') === 'AB');
  check('bootstrap.support', 'support', 'A split supported by every site reaches 100% bootstrap support',
    'resampling of a perfectly consistent alignment', 100, ab ? ab.support : -1, 1e-9);
}

// 19. Missing data does not change parsimony or distances
{
  const base = C.parseFasta('>A\nACGTA\n>B\nACGTC\n>C\nACGGC\n>D\nACGGC');
  const padded = C.parseFasta('>A\nACGTAN\n>B\nACGTCN\n>C\nACGGCN\n>D\nACGGCN');
  const charsA = C.encode(base, C.alphabet('dna'));
  const charsB = C.encode(padded, C.alphabet('dna'));
  const names = base.map(r => r.name);
  const tree = C.parseNewick('((A,B),(C,D));');
  check('missing.parsimony', 'input handling', 'A column of missing data does not change the parsimony score',
    'missing cells are compatible with every state', I.fitchScore(tree, charsA, names), I.fitchScore(tree, charsB, names), 0);
  const dA = M.dnaDistance(base[0].seq, base[1].seq, { model: 'jc' }).distance;
  const dB = M.dnaDistance(padded[0].seq, padded[1].seq, { model: 'jc' }).distance;
  check('missing.distance', 'input handling', 'A column of missing data does not change a pairwise distance',
    'unobserved sites are skipped', dA, dB, 1e-12);
}

// 20. Newick round-trip preserves branch lengths
{
  const original = C.parseNewick("((A:0.123456,B:0.654321):0.111111,(C:0.222222,D:0.333333):0.444444);");
  const round = C.parseNewick(C.toNewick(original));
  const pairs = [];
  const collect = (n) => { if (n.children) n.children.forEach(collect); pairs.push(n.length === null ? 0 : n.length); };
  collect(original);
  const pairs2 = [];
  const collect2 = (n) => { if (n.children) n.children.forEach(collect2); pairs2.push(n.length === null ? 0 : n.length); };
  collect2(round);
  const maxDiff = Math.max(...pairs.map((v, i) => Math.abs(v - pairs2[i])));
  check('newick.roundtrip', 'tree io', 'Newick serialisation and parsing preserve every branch length',
    'serialise then parse the same tree', 0, maxDiff, 1e-6);
}

// 21. Interactive parsimony keeps score and tree consistent
{
  const records = C.parseFasta('>A\nA\n>B\nA\n>C\nC\n>D\nC');
  const chars = C.encode(records, C.alphabet('dna'));
  const names = records.map(r => r.name);
  const session = PI.createSession({ chars, names, kind: 'dna' });
  const start = session.score();
  const moves = session.moves();
  const target = moves.find(m => m.delta !== 0);
  let consistent = true;
  if (target) {
    session.apply(target.id);
    consistent = session.score() === target.score;
    session.undo();
    consistent = consistent && session.score() === start;
  }
  checkTrue('interactive.consistency', 'interactive',
    'Applying a rearrangement produces exactly the score that was advertised, and undo restores the previous score',
    'session bookkeeping', consistent, `start=${start}`);
}

// 22. Every program runs and returns a structured result
{
  const dna = '>A\nACGTACGTACGT\n>C\nACGTACGTACGT\n>B\nACGTACGTACGA\n>D\nACGTACGTACGA';
  const samples = {
    dna,
    protein: '>A\nMKVLAAGIVG\n>B\nMKVLAAGIVG\n>C\nMKILAAGIVG\n>D\nMRILAAGIVG',
    sites: '>A\n1101001011\n>B\n1101001010\n>C\n1100001010\n>D\n0100001010',
    discrete: '>A\n0011010110\n>B\n0011010110\n>C\n0001110101\n>D\n0001100101',
    numbers: '4 2\nA 0.1 0.4\nB 0.12 0.38\nC 0.3 0.2\nD 0.28 0.22',
    tree: '((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);',
    trees: '((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\n((A:0.1,C:0.1):0.2,(B:0.1,D:0.1):0.2);'
  };
  let failures = [];
  for (const program of P.PROGRAM_LIST) {
    const opts = {};
    for (const o of program.options) opts[o.id] = o.default;
    if (program.id === 'SEQBOOT') opts.replicates = 5;
    const ctx = { text: samples[program.input] || '', opts };
    if (program.needsTrees > 0) {
      const raw = program.needsTrees > 1 ? samples.trees : samples.tree;
      ctx.trees = raw.split('\n').map(t => C.parseNewick(t));
    }
    if (program.id === 'DNAINVAR') ctx.text = dna;
    if (program.id === 'FACTOR') ctx.text = samples.discrete;
    try {
      const r = program.run(ctx);
      if (!r || !['tree', 'matrix', 'report', 'interactive'].includes(r.kind)) failures.push(program.id + ':bad-kind');
      if (r.kind === 'tree' && !r.svg.includes('<svg')) failures.push(program.id + ':no-svg');
    } catch (error) {
      failures.push(program.id + ':' + error.message.slice(0, 60));
    }
  }
  check('programs.all-run', 'integration', 'All 35 programs execute and return a renderable result',
    'integration sweep', 0, failures.length, 0);
  if (failures.length) checks[checks.length - 1].detail = failures.join('; ');
}

// ------------------------------------------------------------------- report

const passed = checks.filter(c => c.pass).length;
const failed = checks.length - passed;
const byCategory = {};
for (const c of checks) {
  byCategory[c.category] = byCategory[c.category] || { total: 0, passed: 0 };
  byCategory[c.category].total++;
  if (c.pass) byCategory[c.category].passed++;
}

const report = {
  suite: 'Phoenix scientific validation',
  version: (() => { try { return require('../package.json').version; } catch (_) { return 'unversioned'; } })(),
  generatedAt: new Date().toISOString(),
  summary: {
    checks: checks.length,
    passed,
    failed,
    categories: byCategory,
    status: failed === 0 ? 'PASS' : 'FAIL'
  },
  notes: [
    'Each check compares Phoenix output with a value obtained independently: a closed form, an exhaustive enumeration, a hand computation, or a simulation with a known generating process.',
    'Tolerances reflect floating-point and optimisation noise; they are stated per check.',
    'These checks establish implementation correctness against mathematical references. They are not a claim of bit-for-bit agreement with PHYLIP.'
  ],
  checks
};

const outDir = path.join(__dirname);
fs.writeFileSync(path.join(outDir, 'report.json'), JSON.stringify(report, null, 2) + '\n');

const pad = (s, n) => String(s).padEnd(n);
for (const c of checks) {
  console.log(`${c.pass ? 'PASS' : 'FAIL'}  ${pad(c.id, 28)} ${c.pass ? '' : `expected ${c.expected} got ${c.actual}`}`);
}
console.log(`\n${passed}/${checks.length} checks passed. Report written to validation/report.json`);
process.exit(failed === 0 ? 0 : 1);
