#!/usr/bin/env node
/* validation/benchmark.js — timing of representative analyses.
 * Run: node validation/benchmark.js   (writes validation/benchmark.json)
 */
'use strict';
const fs = require('fs');
const os = require('os');
const C = require('../phoenix-core.js');
const M = require('../phoenix-models.js');
const I = require('../phoenix-infer.js');
const S = require('../phoenix-support.js');
const P = require('../phoenix-programs.js');

function time(fn) {
  const t0 = process.hrtime.bigint();
  const value = fn();
  const t1 = process.hrtime.bigint();
  return { ms: Number(t1 - t0) / 1e6, value };
}

function simulatedAlignment(nTaxa, nSites, seed) {
  const rand = C.rng(seed);
  const model = M.buildQ({ kind: 'dna', model: 'jc' });
  // Star-ish tree with one internal node; enough structure to be realistic.
  const seqs = Array.from({ length: nTaxa }, () => []);
  for (let s = 0; s < nSites; s++) {
    let root = Math.floor(rand() * 4);
    for (let t = 0; t < nTaxa; t++) {
      const P_ = M.transitionMatrix(model.Q, 0.08, 4);
      let r = rand(), st = 3;
      for (let j = 0; j < 4; j++) { r -= P_[root][j]; if (r <= 0) { st = j; break; } }
      seqs[t].push('ACGT'[st]);
    }
    void root;
    root = 0;
  }
  return seqs.map((s, i) => ({ name: `T${i + 1}`, seq: s.join('') }));
}

const rows = [];
function bench(name, detail, fn) {
  const r = time(fn);
  rows.push({ name, detail, ms: Number(r.ms.toFixed(1)) });
  console.log(`${name.padEnd(34)} ${String(r.ms.toFixed(1)).padStart(9)} ms   ${detail}`);
}

for (const [nTaxa, nSites] of [[10, 500], [30, 1000]]) {
  const records = simulatedAlignment(nTaxa, nSites, 12345);
  const fasta = records.map(r => `>${r.name}\n${r.seq}`).join('\n');
  bench('DNADIST (JC)', `${nTaxa} taxa × ${nSites} sites`, () => P.get('DNADIST').run({ text: fasta, opts: { model: 'jc' } }));
  bench('NEIGHBOR (NJ)', `${nTaxa} taxa from the same matrix`, () => P.get('NEIGHBOR').run({ text: fasta, opts: { method: 'nj', model: 'jc' } }));
  bench('DNAPARS', `${nTaxa} taxa × ${nSites} sites`, () => P.get('DNAPARS').run({ text: fasta, opts: { restarts: 2 } }));
}

for (const nTaxa of [8, 16]) {
  const records = simulatedAlignment(nTaxa, 300, 999);
  const fasta = records.map(r => `>${r.name}\n${r.seq}`).join('\n');
  bench('DNAML (HKY)', `${nTaxa} taxa × 300 sites`, () => P.get('DNAML').run({ text: fasta, opts: { model: 'hky', kappa: 2 } }));
}

{
  const records = simulatedAlignment(12, 400, 777);
  const fasta = records.map(r => `>${r.name}\n${r.seq}`).join('\n');
  bench('SEQBOOT + parsimony + CONSENSE', '12 taxa × 400 sites, 50 replicates', () =>
    P.get('SEQBOOT').run({ text: fasta, opts: { replicates: 50, method: 'bootstrap', seed: 5, summarise: true, treeMethod: 'parsimony', rule: 'majority' } }));
}

const report = {
  suite: 'Phoenix performance benchmark',
  generatedAt: new Date().toISOString(),
  environment: {
    runtime: `node ${process.version}`,
    platform: `${os.platform()} ${os.arch()}`,
    cpu: (os.cpus()[0] || {}).model || 'unknown',
    cores: os.cpus().length,
    note: 'Containerised CI sandbox; absolute timings will differ on other hardware.'
  },
  results: rows
};
fs.writeFileSync(require('path').join(__dirname, 'benchmark.json'), JSON.stringify(report, null, 2) + '\n');
console.log('\nWrote validation/benchmark.json');
