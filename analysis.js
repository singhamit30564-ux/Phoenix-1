/* PHOENIX's small, transparent distance-based phylogeny engine. */
(function (root, factory) {
  const api = factory();
  if (typeof module === 'object' && module.exports) module.exports = api;
  else root.Phylo = api;
})(typeof globalThis !== 'undefined' ? globalThis : this, function () {
  function parseFasta(text) {
    const records = []; let current = null;
    for (const raw of text.split(/\r?\n/)) {
      const line = raw.trim(); if (!line) continue;
      if (line.startsWith('>')) {
        const name = line.slice(1).trim();
        if (!name) throw new Error('Every FASTA header needs a name.');
        current = { name, sequence: '' }; records.push(current);
      } else {
        if (!current) throw new Error('Sequences must start with a FASTA header (>name).');
        current.sequence += line.replace(/\s/g, '').toUpperCase();
      }
    }
    if (records.length < 2) throw new Error('Please provide at least two sequences.');
    if (records.length > 100) throw new Error('Please provide no more than 100 sequences.');
    if (new Set(records.map(r => r.name)).size !== records.length) throw new Error('Sequence names must be unique.');
    if (!records[0].sequence.length) throw new Error('Sequences cannot be empty.');
    if (records.some(r => r.sequence.length !== records[0].sequence.length)) throw new Error('Sequences must be aligned and have equal lengths.');
    if (records.some(r => /[^ACGTN?-]/.test(r.sequence))) throw new Error('Only DNA bases A, C, G, T and missing characters N, ? or - are supported.');
    return records;
  }
  function distance(a, b, model) {
    let compared = 0, differences = 0;
    for (let i = 0; i < a.length; i++) {
      if (!'ACGT'.includes(a[i]) || !'ACGT'.includes(b[i])) continue;
      compared++; if (a[i] !== b[i]) differences++;
    }
    if (!compared) throw new Error('A pair of sequences has no comparable DNA sites.');
    const p = differences / compared;
    if (model === 'jc' && p >= 0.75) throw new Error('Jukes–Cantor distance is undefined when 75% or more of comparable sites differ.');
    return model === 'jc' ? -0.75 * Math.log(1 - 4 * p / 3) : p;
  }
  function analyze(text, model = 'p') {
    if (!['p', 'jc'].includes(model)) throw new Error('Unknown distance model.');
    const records = parseFasta(text);
    const matrix = records.map(a => records.map(b => a === b ? 0 : distance(a.sequence, b.sequence, model)));
    // UPGMA joins the closest clusters using size-weighted average distances.
    let clusters = records.map((r, i) => ({ name: r.name, size: 1, height: 0, index: i }));
    let d = matrix.map(row => row.slice());
    while (clusters.length > 1) {
      let x = 0, y = 1;
      for (let i = 0; i < clusters.length; i++) for (let j = i + 1; j < clusters.length; j++) if (d[i][j] < d[x][y]) { x = i; y = j; }
      const a = clusters[x], b = clusters[y], height = d[x][y] / 2;
      const merged = { children: [a, b], height, size: a.size + b.size };
      const rest = clusters.map((_, i) => i).filter(i => i !== x && i !== y);
      const nextD = rest.map(i => rest.map(j => d[i][j]));
      const averages = rest.map(i => (d[x][i] * a.size + d[y][i] * b.size) / merged.size);
      nextD.forEach((row, i) => row.push(averages[i])); nextD.push([...averages, 0]);
      clusters = [...rest.map(i => clusters[i]), merged]; d = nextD;
    }
    return { records, matrix, tree: clusters[0], sites: records[0].sequence.length, model };
  }
  function newick(node, parentHeight) {
    const label = node.children ? '(' + node.children.map(c => newick(c, node.height)).join(',') + ')' : "'" + node.name.replace(/'/g, "''") + "'";
    return label + (parentHeight === undefined ? ';' : ':' + Math.max(0, parentHeight - node.height).toFixed(6));
  }
  return { parseFasta, distance, analyze, newick };
});
