/* demos.js — three guided end-to-end workflows.
 *
 * Each demo runs real analyses through the same programs the suite exposes;
 * nothing here fabricates results. Step outputs are collected and rendered by
 * the program-suite output pane.
 */
(function () {
  'use strict';

  const C = window.PCore;
  const P = window.PPrograms;
  const PI = window.PInteractive;
  const S = window.PSupport;

  // Synthetic datasets written for Phoenix. No third-party data is bundled.
  const DATA = {
    primates: `>Human
ATGCTAGCTAGCTACGATCGATCGTACGATCG
>Chimpanzee
ATGCTAGCTAGCTACGATCGATCGTACGATCA
>Gorilla
ATGCTAGCTAGCTACGATCGATCATACGATCA
>Orangutan
ATGCTAGCTAGCTACGATCAATCATACGATCA
>Gibbon
ATGCTAGCTAGTTACGATCAATCATACGATCA`,
    characters: `>TaxonA
001101011001
>TaxonB
001101011001
>TaxonC
000111010101
>TaxonD
000110010102`,
    beetles: `>BeetleA
TTAAGACTCAGTGTCCACGACAAGCGCTGTCACACTTTCACTGAAACTACGGACGACACCCCGGCGCTGAGTGTAGGTAGGGCTAACCACGGTAGATAGAGGTGTACCAGTCGAATGAATCAGCTCACTTCGGACGACGAAAGCTAGTTAATGCTCCTACGGCGGTCGTAGCCAAAGCCAAACTACATCTGTCTGGGTACGACATTCCATACAAATGCCA
>BeetleB
TTCAGACTCAGTGCCCACGACAAGCGCTGTCACACTTTCACTGAAACTACGGACCACACCCCCGCGCTGAGTGTAGGTTGGACTAACCATGGTCGATAGAGATGTACCTGTCGAATGACTCAGCTCACTTCGGACGACGCAAGCTAGTTACTGCTCCTACGGCGGTTGAAGCCAAAGCCCAACTACCTCTGTCTCGGTACTACATTCCATACAAAGGCCA
>BeetleC
TTCATTCTCAGCGTAATTTACAAGCTCTGTCTCTCCAGCTCTGAAACTTCGGAGTACAAACCGGCACTGTGTAAAGGCAGGTCTTAGCATGGTAGATATAGGTGTACTGGTCGAATGAATTGGCTCACCTCCGCAGACTGAACGTTGAGAACGCGCCAACGGTGATGAATGCCAAAGTCCAACTAAATCTGTCTCCGTACAACATTGGATACAAAGGCAA
>BeetleD
TTCATTCTCAGCGTTCTGTACAAGCTCTGTCTCTCCATCACTGAAACTACGGTGTACAAACCGGCACTGTGTAAAGGCAGGGCTGAATATAGTAGATATAGGTGTTCAGGTCGAATGAATGAGCTCACCTCCGAAGACTGAACATAGTGAAGGCGCCTACGGTGATGAATACCAAAGTCCAATAAAATCTGTCTCCGTACAACATTGGATACAAAGGCAA
>BeetleE
TTCCATCTACAGGTTCTCTAAAAGCGAAGACACATTTTCCCTGGATATGCACACGTCAACCCGCGACTACATTTAGGCAGGGCGAGAAATGGGCAACTGAGATGTTCATAAGGTCTTAAGGACTTCACCTGAGACGACTAAAGCTTAAGAATGCTCCTAGGGAGTCTAATGGCAGAACTGATCTAAGACTGACTGGGTACGTGATTCTAGACACAGGCAA
>BeetleF
TTCCATCTACAGCTTCTCTAGAAGCGATGACACACTTTCCCTGGATATGGACACGTCAACCCGCGACTACATTTAGGCAGGGCGATAAAAGGGCAACTGAGATGTTCATAAGGTCTTAAGGAATTCACCTGAGAAGACTTAAGCTTAAGAATGCTCCTAGGGAGTCTAATGGCAGAACTGATTCAAGACTGACTGGGTACGTGAATATAGACACAGGCAA`
  };

  function runProgram(id, text, opts) {
    const program = P.get(id);
    const options = {};
    for (const o of program.options) options[o.id] = o.default;
    Object.assign(options, opts || {});
    return program.run({ text, opts: options });
  }

  function step(title, programId, result, caption) {
    return { title, programId, result, caption };
  }

  const DEMOS = [
    {
      id: 'dna',
      label: 'A',
      title: 'DNA distances and a tree',
      blurb: 'Turn an aligned DNA alignment into a Jukes–Cantor distance matrix and a neighbour-joining tree.',
      dataset: 'primates',
      dataLabel: '5 primate sequences, 32 sites (synthetic)',
      minutes: 'seconds',
      build() {
        const distance = runProgram('DNADIST', DATA.primates, { model: 'jc' });
        const tree = runProgram('NEIGHBOR', DATA.primates, { method: 'nj', model: 'jc' });
        return {
          title: 'DNA distances and a tree',
          blurb: 'Jukes–Cantor distances between five primate sequences, then neighbour joining.',
          steps: [
            step('1. Pairwise distances', 'DNADIST', distance,
              'The Jukes–Cantor correction converts the observed proportion of differing sites into substitutions per site.'),
            step('2. Tree from the distances', 'NEIGHBOR', tree,
              'Neighbour joining builds a tree from that matrix. Branch lengths are in substitutions per site.')
          ],
          tree: tree.tree,
          text: `${distance.text}\n\n${tree.ascii}`,
          notes: [
            'Distances are computed from the alignment, not from the tree, so the two steps can be checked against each other.',
            'Neighbour joining is unrooted; Phoenix displays it rooted on the first internal node.'
          ]
        };
      }
    },
    {
      id: 'parsimony',
      label: 'B',
      title: 'Parsimony and ancestral states',
      blurb: 'Score a discrete-character matrix with Wagner parsimony and read the implied state changes.',
      dataset: 'characters',
      dataLabel: '4 taxa, 12 discrete characters (synthetic)',
      minutes: 'seconds',
      build() {
        const parsed = runProgram('PARS', DATA.characters, { restarts: 6 });
        const records = C.parseInput(DATA.characters);
        const chars = C.encode(records, C.alphabet('discrete'));
        const names = records.map(r => r.name);
        const session = PI.createSession({
          chars, names, kind: 'discrete', weights: { wagner: 1 }, initialTree: parsed.tree
        });
        const rows = session.siteReconstruction(12);
        const table = [
          'Character   cost   branches changed   root state',
          ...rows.map(r => `${String(r.site).padStart(9)}${String(r.cost).padStart(7)}${String(r.changes).padStart(19)}   ${r.rootStates.join(',') || '?'}`)
        ].join('\n');
        return {
          title: 'Parsimony and ancestral states',
          blurb: 'Wagner parsimony on ordered multistate characters, with the most-parsimonious reconstruction of each character.',
          steps: [
            step('1. Wagner parsimony', 'PARS', parsed,
              'States are ordered, so a change of k states costs k steps. The score is the total number of steps.'),
            step('2. Ancestral reconstruction', 'MOVE', {
              kind: 'report', text: table, stats: [], notes: []
            }, 'The reconstruction uses the same criterion as the score, so the per-character costs sum to it.')
          ],
          tree: parsed.tree,
          text: `${parsed.ascii}\n\n${table}`,
          notes: [
            'Wagner parsimony treats states as ordered; unordered data should use DNAPARS or PROTPARS instead.',
            'Equally parsimonious reconstructions can exist; ties are broken toward the lowest state index.'
          ],
          openWith: { program: 'MOVE', text: DATA.characters }
        };
      }
    },
    {
      id: 'bootstrap',
      label: 'C',
      title: 'Bootstrap support',
      blurb: 'Resample a simulated alignment, rebuild trees and summarise them as a consensus. The generating tree is known, so the support values can be checked against the truth.',
      dataset: 'beetles',
      dataLabel: '6 beetle sequences, 220 sites, simulated under JC69',
      minutes: 'a few seconds',
      build() {
        const boot = runProgram('SEQBOOT', DATA.beetles, {
          replicates: 50, method: 'bootstrap', seed: 20240301,
          summarise: true, treeMethod: 'parsimony', rule: 'majority'
        });
        return {
          title: 'Bootstrap support',
          blurb: 'Fifty bootstrap replicates searched with parsimony, summarised by majority-rule consensus.',
          steps: [
            step('1. Resample, search and summarise', 'SEQBOOT', boot,
              'Columns are resampled with replacement, each replicate is searched, and splits present in most replicate trees are kept.')
          ],
          tree: boot.tree,
          text: boot.ascii,
          notes: [
            'Support values are percentages of replicate trees containing that split. They describe the data, not a posterior probability.',
            'The seed is fixed, so the same run reproduces exactly; change the seed in the program options to see sampling variation.'
          ]
        };
      }
    }
  ];

  function run(demoId) {
    const demo = DEMOS.find(d => d.id === demoId);
    if (!demo) throw new Error(`Unknown demo "${demoId}".`);
    const built = demo.build();
    return Object.assign({
      kind: 'demo',
      demoId: demo.id,
      dataset: demo.dataLabel,
      data: DATA[demo.dataset]
    }, built);
  }

  function mount() {
    const host = document.getElementById('demoGrid');
    if (!host) return;
    host.innerHTML = DEMOS.map(d => `
      <article class="demo-card">
        <div class="demo-head">
          <span class="demo-letter">${d.label}</span>
          <h3>${d.title}</h3>
        </div>
        <p>${d.blurb}</p>
        <dl class="demo-meta">
          <div><dt>DATA</dt><dd>${d.dataLabel}</dd></div>
          <div><dt>RUNS</dt><dd>${d.minutes}</dd></div>
        </dl>
        <button class="demo-run" data-demo="${d.id}">Run demo <span>↗</span></button>
      </article>`).join('');
    host.querySelectorAll('[data-demo]').forEach(btn => {
      btn.addEventListener('click', () => {
        const id = btn.dataset.demo;
        const label = btn.innerHTML;
        btn.disabled = true;
        btn.innerHTML = 'Running…';
        window.setTimeout(() => {
          try {
            const result = run(id);
            window.PhoenixUI.showResult(result, 'demo');
            document.getElementById('suite-output').scrollIntoView({ behavior: 'smooth', block: 'start' });
          } catch (error) {
            window.PhoenixUI.showError(`Demo ${id} failed: ${error.message}`);
          } finally {
            btn.disabled = false;
            btn.innerHTML = label;
          }
        }, 20);
      });
    });
  }

  window.PhoenixDemos = { DEMOS, run, mount, DATA };
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', mount);
  else mount();
})();
