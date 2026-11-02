const { test } = require('node:test');
const assert = require('node:assert/strict');
const { parseFasta, distance, analyze, newick } = require('../analysis.js');
test('FASTA parsing and validation', () => {
  assert.deepEqual(parseFasta('>a\nac gt\n>b\nACGT').map(x => x.sequence), ['ACGT','ACGT']);
  assert.throws(() => parseFasta('>a\nAC\n>b\nA'), /equal lengths/);
  assert.throws(() => parseFasta('>a\nAC\n>a\nAC'), /unique/);
  assert.throws(() => parseFasta('>a\nAX\n>b\nAC'), /Only DNA/);
});
test('pairwise distance and missing sites', () => {
  assert.equal(distance('ACN-', 'ATGG', 'p'), .5);
  assert.throws(() => distance('NN', '--', 'p'), /no comparable/);
  assert.throws(() => distance('AAAA', 'TTTA', 'jc'), /undefined/);
  assert.ok(Math.abs(distance('AAAA', 'AAAT', 'jc') - 0.304098831) < 1e-8);
});
test('UPGMA clusters closest sequences and exports escaped names', () => {
  const r = analyze(">a's\nAAAA\n>b\nAAAT\n>c\nTTTT");
  assert.equal(r.sites, 4);
  assert.equal(newick(r.tree), "('c':0.437500,('a''s':0.125000,'b':0.125000):0.312500);");
});
