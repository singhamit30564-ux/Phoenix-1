#!/usr/bin/env node
/* validation/generate-compatibility.js
 * Writes docs/COMPATIBILITY.md from the live program registry, so the
 * compatibility matrix can never drift away from the code that implements it.
 */
'use strict';

const fs = require('fs');
const path = require('path');
const P = require('../phoenix-programs.js');

const ORDER = [
  ['Parsimony and compatibility', ['PROTPARS', 'DNAPARS', 'DNAPENNY', 'DNAMOVE', 'DNACOMP', 'PARS', 'MIX', 'PENNY', 'MOVE', 'DOLLOP', 'DOLPENNY', 'DOLMOVE', 'CLIQUE']],
  ['Likelihood', ['DNAML', 'DNAMLK', 'PROML', 'PROMLK', 'RESTML']],
  ['Distances', ['DNADIST', 'PROTDIST', 'RESTDIST', 'GENDIST', 'TREEDIST']],
  ['Distance trees', ['NEIGHBOR', 'FITCH', 'KITSCH', 'CONTML', 'CONTRAST']],
  ['Resampling and diagnostics', ['SEQBOOT', 'CONSENSE', 'RETREE', 'DNAINVAR', 'FACTOR']],
  ['Drawing', ['DRAWGRAM', 'DRAWTREE']]
];

const LEGEND = {
  validated: 'Fully implemented, with numerical checks against an independent reference.',
  'documented differences': 'Correct method, but a feature or optimisation of the original differs and is documented.',
  partial: 'Part of the original behaviour is missing.',
  exploratory: 'A reasoned approach whose numbers should not be treated as reference values.'
};

const rows = [];
for (const [group, ids] of ORDER) {
  rows.push(`\n### ${group}\n`);
  rows.push('| Program | Status | What is implemented | Difference from PHYLIP |');
  rows.push('| --- | --- | --- | --- |');
  for (const id of ids) {
    const program = P.get(id);
    const meta = P.statusOf(id);
    rows.push(`| \`${id}\` | ${meta.status} | ${program.summary.replace(/\|/g, '\\|')} | ${meta.ref.replace(/\|/g, '\\|')} |`);
  }
}

const counts = P.catalogue().reduce((acc, p) => {
  acc[p.status] = (acc[p.status] || 0) + 1;
  return acc;
}, {});

const doc = `# PHYLIP compatibility matrix

This file is generated from the program registry by
\`node validation/generate-compatibility.js\`. Do not edit it by hand; edit
\`phoenix-programs.js\` and regenerate.

Phoenix is an independent reimplementation. It is **not** a port of the PHYLIP
source, and no program here claims bit-for-bit agreement with PHYLIP output.
Status labels mean:

${Object.entries(LEGEND).map(([k, v]) => `- **${k}** — ${v}`).join('\n')}

## Summary

| Status | Programs |
| --- | --- |
${Object.entries(LEGEND).map(([k]) => `| ${k} | ${counts[k] || 0} |`).join('\n')}
| **total** | **${P.PROGRAM_LIST.length}** |

Only a program whose method is checked numerically by
\`validation/validate.js\` is labelled **validated**.
${rows.join('\n')}
`;

const target = path.join(__dirname, '..', 'docs', 'COMPATIBILITY.md');
fs.mkdirSync(path.dirname(target), { recursive: true });
fs.writeFileSync(target, doc);
console.log(`Wrote ${path.relative(process.cwd(), target)} (${P.PROGRAM_LIST.length} programs)`);
