/* programs.js — the program browser: catalogue, option forms, run and render. */
(function () {
  'use strict';
  const $ = id => document.getElementById(id);
  const P = window.PPrograms;
  const C = window.PCore;
  const S = window.PSupport;

  const EXAMPLES = {
    dna: `>Human
ATGCTAGCTAGCTACGATCGATCGTACGATCG
>Chimpanzee
ATGCTAGCTAGCTACGATCGATCGTACGATCA
>Gorilla
ATGCTAGCTAGCTACGATCGATCATACGATCA
>Orangutan
ATGCTAGCTAGCTACGATCAATCATACGATCA
>Gibbon
ATGCTAGCTAGTTACGATCAATCATACGATCA`,
    protein: `>Alpha
MKVLAAGIVGLNLGGKST
>Beta
MKVLAAGIVGLNLGGKST
>Gamma
MKILAAGIVGLNIGGKST
>Delta
MRILAAGIVGLNIGGKTT`,
    sites: `>Sample1
11010010110100101101
>Sample2
11010010110100101100
>Sample3
11000010110100101100
>Sample4
01000010110000101100`,
    discrete: `>Taxon1
001101011001
>Taxon2
001101011001
>Taxon3
000111010101
>Taxon4
000110010102`,
    numbers: `5 4
A 0.10 0.40 0.20 0.30
B 0.12 0.38 0.22 0.28
C 0.30 0.20 0.25 0.25
D 0.28 0.22 0.26 0.24
E 0.50 0.10 0.20 0.20`,
    tree: `((Human:0.08,Chimpanzee:0.09):0.05,(Gorilla:0.14,Orangutan:0.16):0.07,Gibbon:0.21);`,
    trees: `((Human:0.08,Chimpanzee:0.09):0.05,(Gorilla:0.14,Orangutan:0.16):0.07,Gibbon:0.21);
((Human:0.07,Gorilla:0.09):0.06,(Chimpanzee:0.15,Orangutan:0.17):0.05,Gibbon:0.22);
((Human:0.07,Orangutan:0.09):0.06,(Gorilla:0.15,Chimpanzee:0.17):0.05,Gibbon:0.22);
((Human:0.09,Chimpanzee:0.07):0.04,(Gorilla:0.13,Orangutan:0.15):0.08,Gibbon:0.20);`
  };

  const INPUT_LABEL = {
    dna: { title: 'DNA sequences', hint: 'FASTA or PHYLIP · A C G T N ? -', kind: 'dna' },
    protein: { title: 'Protein sequences', hint: 'FASTA or PHYLIP · 20 amino acids', kind: 'protein' },
    sites: { title: 'Restriction sites', hint: 'FASTA or PHYLIP · 1 = site present, 0 = absent', kind: 'sites' },
    discrete: { title: 'Discrete characters', hint: 'FASTA or PHYLIP · states 0–9', kind: 'discrete' },
    numbers: { title: 'Numeric matrix', hint: 'name followed by numbers, one taxon per line', kind: 'numbers' },
    tree: { title: 'Tree', hint: 'Newick', kind: 'tree' },
    trees: { title: 'Trees', hint: 'One Newick tree per line', kind: 'tree' }
  };

  const store = Object.assign({}, EXAMPLES);
  let current = null;
  let last = null;
  let activeInput = null;

  // ------------------------------------------------------------- catalogue

  function ensureStatusNodes() {
    if ($('programStatus')) return;
    const badge = document.createElement('span');
    badge.id = 'programStatus';
    badge.className = 'status-badge';
    badge.hidden = true;
    $('programCategory').insertAdjacentElement('afterend', badge);
    const ref = document.createElement('p');
    ref.id = 'programReference';
    ref.className = 'program-reference';
    $('programCategory').closest('.card-heading').insertAdjacentElement('afterend', ref);
  }

  function renderCatalogue(filter) {
    const container = $('programList');
    const term = (filter || '').trim().toLowerCase();
    const groups = new Map();
    for (const p of P.catalogue()) {
      if (term && !(`${p.id} ${p.title} ${p.category} ${p.summary}`.toLowerCase().includes(term))) continue;
      if (!groups.has(p.category)) groups.set(p.category, []);
      groups.get(p.category).push(p);
    }
    if (!groups.size) {
      container.innerHTML = '<p class="suite-empty">No program matches that search.</p>';
      return;
    }
    const parts = [];
    for (const category of P.CATEGORIES) {
      if (!groups.has(category)) continue;
      parts.push(`<div class="program-group"><span>${category.toUpperCase()}</span>`);
      for (const p of groups.get(category)) {
        parts.push(`<button class="program-item${current && current.id === p.id ? ' is-active' : ''}" data-program="${p.id}" title="${escapeAttr(p.summary)}">
          <b>${p.id}</b><i>${escapeHtml(p.title)}</i></button>`);
      }
      parts.push('</div>');
    }
    container.innerHTML = parts.join('');
    container.querySelectorAll('[data-program]').forEach(btn => {
      btn.addEventListener('click', () => select(btn.dataset.program));
    });
  }

  function escapeHtml(s) {
    return String(s).replace(/[<>&]/g, c => ({ '<': '&lt;', '>': '&gt;', '&': '&amp;' }[c]));
  }
  function escapeAttr(s) {
    return escapeHtml(s).replace(/"/g, '&quot;');
  }

  // -------------------------------------------------------------- selection

  function select(id) {
    const program = P.PROGRAM_LIST.find(p => p.id === id);
    if (!program) return;
    current = program;
    last = null;
    renderCatalogue($('programSearch').value);
    $('programName').textContent = program.id;
    $('programTitle').textContent = program.title;
    $('programCategory').textContent = program.category.toUpperCase();
    const meta = P.statusOf(program.id);
    const badge = $('programStatus');
    if (meta) {
      badge.textContent = meta.status.toUpperCase();
      badge.className = 'status-badge status-' + meta.status.replace(/\s+/g, '-');
      badge.title = meta.ref || '';
      badge.hidden = false;
    } else {
      badge.hidden = true;
    }
    $('programReference').textContent = meta && meta.ref ? meta.ref : '';
    $('programSummary').textContent = program.summary;
    renderOptions(program);
    configureInput(program);
    clearOutput('Run ' + program.id + ' to see results.');
    $('programError').textContent = '';
  }

  function renderOptions(program) {
    const box = $('programOptions');
    if (!program.options.length) {
      box.innerHTML = '<p class="no-options">This program takes no options.</p>';
      return;
    }
    box.innerHTML = program.options.map(opt => {
      const id = 'opt_' + opt.id;
      if (opt.type === 'select') {
        return `<div><label for="${id}">${opt.label}</label>
          <select id="${id}" data-option="${opt.id}">
            ${opt.choices.map(c => `<option value="${c}"${c === opt.default ? ' selected' : ''}>${c}</option>`).join('')}
          </select></div>`;
      }
      if (opt.type === 'checkbox') {
        return `<div class="check-option"><label for="${id}"><input id="${id}" type="checkbox" data-option="${opt.id}"${opt.default ? ' checked' : ''}/> ${opt.label}</label></div>`;
      }
      if (opt.type === 'text') {
        return `<div><label for="${id}">${opt.label}</label>
          <input id="${id}" type="text" data-option="${opt.id}" value="${escapeAttr(opt.default || '')}" /></div>`;
      }
      return `<div><label for="${id}">${opt.label}</label>
        <input id="${id}" type="number" data-option="${opt.id}" value="${opt.default}" min="${opt.min}" max="${opt.max}" step="${opt.step || 1}" /></div>`;
    }).join('');
  }

  function configureInput(program) {
    const treeOnly = program.input === 'tree' || program.input === 'trees';
    const info = INPUT_LABEL[program.input] || INPUT_LABEL.dna;
    $('programInputWrap').hidden = treeOnly;
    $('programTreeWrap').hidden = !(program.needsTrees > 0);
    if (treeOnly) $('programTreeWrap').hidden = false;
    $('programInputLabel').textContent = info.title;
    $('programInputHint').textContent = info.hint;
    if (!treeOnly) {
      activeInput = program.input;
      if (!store[program.input]) store[program.input] = EXAMPLES.dna;
      $('programInput').value = store[program.input];
    }
    if (program.needsTrees > 0) {
      const key = program.needsTrees > 1 ? 'trees' : 'tree';
      if (!store[key]) store[key] = EXAMPLES[key] || EXAMPLES.tree;
      $('programTree').value = store[key];
      $('programTreeLabel').textContent = program.needsTrees > 1
        ? `Trees — ${program.needsTrees} or more, one Newick tree per line`
        : 'Tree — Newick';
    }
  }

  function readOptions() {
    const opts = {};
    $('programOptions').querySelectorAll('[data-option]').forEach(el => {
      const key = el.dataset.option;
      if (el.type === 'checkbox') opts[key] = el.checked;
      else if (el.type === 'number') opts[key] = el.value === '' ? undefined : Number(el.value);
      else opts[key] = el.value;
    });
    return opts;
  }

  // ------------------------------------------------------------------- run

  function run() {
    if (!current) return;
    const button = $('runProgram');
    button.disabled = true;
    $('programError').textContent = '';
    const started = performance.now();
    requestAnimationFrame(() => {
      try {
        if (activeInput !== null && $('programInput')) store[activeInput] = $('programInput').value;
        const opts = readOptions();
        const ctx = { text: activeInput !== null ? $('programInput').value : '', opts };
        if (current.needsTrees > 0) {
          const raw = $('programTree').value;
          store[current.needsTrees > 1 ? 'trees' : 'tree'] = raw;
          ctx.trees = raw.split(/;\s*(?:\r?\n|$)/).map(s => s.trim()).filter(Boolean).map(s => C.parseNewick(s + ';'));
          if (ctx.trees.length < (current.needsTrees || 1)) {
            throw new Error(`Supply at least ${current.needsTrees} tree${current.needsTrees > 1 ? 's' : ''}, separated by semicolons.`);
          }
        }
        const result = current.run(ctx);
        result.elapsed = Math.round(performance.now() - started);
        last = result;
        render(result);
      } catch (error) {
        last = null;
        clearOutput('The run stopped before producing a result.');
        $('programError').textContent = error.message;
      } finally {
        button.disabled = false;
      }
    });
  }

  function clearOutput(message) {
    $('programView').hidden = false;
    $('programView').classList.remove('is-tree');
    $('programView').innerHTML = `<div class="output-empty"><span>✳</span>${escapeHtml(message)}</div>`;
    $('programText').textContent = '';
    $('programStats').innerHTML = '';
    $('programNotes').innerHTML = '';
    $('downloadNewick').disabled = true;
    $('downloadText').disabled = true;
  }

  let activeSession = null;

  function render(result) {
    const view = $('programView');
    if (result.kind === 'interactive') {
      activeSession = result.session;
      view.classList.remove('is-tree');
      drawInteractive(result);
      $('downloadNewick').disabled = false;
      $('programNotes').innerHTML = (result.notes || []).map(n => `<li>${escapeHtml(n)}</li>`).join('');
      $('programError').textContent = '';
      return;
    }
    activeSession = null;
    if (result.kind === 'tree') {
      view.innerHTML = result.svg;
      view.classList.add('is-tree');
      view.hidden = false;
    } else {
      view.classList.remove('is-tree');
      view.innerHTML = '';
      view.hidden = true;
    }
    $('programText').textContent = result.kind === 'tree' ? result.ascii : (result.text || '');
    $('programStats').innerHTML = (result.stats || []).map(s =>
      `<div><span>${escapeHtml(s.label)}</span><strong>${escapeHtml(String(s.value))}</strong></div>`).join('')
      + `<div><span>RUNTIME</span><strong>${result.elapsed} ms</strong></div>`;
    $('programNotes').innerHTML = (result.notes || []).map(n => `<li>${escapeHtml(n)}</li>`).join('');
    $('downloadNewick').disabled = result.kind !== 'tree';
    $('downloadText').disabled = !($('programText').textContent || '').trim();
    $('programError').textContent = '';
  }

  function drawInteractive(result) {
    const session = activeSession;
    const score = session.score();
    const moves = session.moves().sort((a, b) => a.delta - b.delta).slice(0, 8);
    const svg = S.drawRectangular(session.tree, { style: 'cladogram', width: 720, rowHeight: 24 });
    const moveButtons = moves.map(m => `
      <button class="move-button${m.delta < 0 ? ' is-improving' : ''}" data-move="${escapeAttr(m.id)}">
        <span>${escapeHtml(m.label)}</span>
        <b>${m.delta === 0 ? '±0' : (m.delta > 0 ? '+' : '') + m.delta}</b>
      </button>`).join('') || '<p class="suite-empty">No rearrangements available.</p>';
    $('programView').classList.remove('is-tree');
    $('programView').innerHTML = `
      <div class="interactive-head">
        <div><span class="mini-label">STEPS ON THIS TREE</span><strong class="interactive-score">${score}</strong></div>
        <div class="interactive-actions">
          <button class="download" id="interactiveUndo" ${session.canUndo() ? '' : 'disabled'}>↶ Undo</button>
          <button class="download" id="interactiveGreedy">↧ Finish search</button>
          <button class="download" id="interactiveReset">↺ New tree</button>
        </div>
      </div>
      <div class="interactive-tree">${svg}</div>
      <div class="interactive-moves"><span class="mini-label">REARRANGEMENTS — APPLY ONE</span>${moveButtons}</div>`;
    $('programText').textContent = siteTable(session);
    $('programStats').innerHTML = (result.stats || []).slice(0, 1).map(st =>
      `<div><span>${escapeHtml(st.label)}</span><strong>${escapeHtml(String(score))}</strong></div>`).join('')
      + `<div><span>CHANGES SINCE START</span><strong>${session.historyLength()}</strong></div>`
      + `<div><span>RUNTIME</span><strong>${result.elapsed} ms</strong></div>`;
    $('downloadText').disabled = false;
    $('programView').querySelectorAll('[data-move]').forEach(btn => {
      btn.addEventListener('click', () => {
        session.apply(btn.dataset.move);
        drawInteractive(result);
      });
    });
    $('interactiveUndo').addEventListener('click', () => { session.undo(); drawInteractive(result); });
    $('interactiveGreedy').addEventListener('click', () => { session.greedy(); drawInteractive(result); });
    $('interactiveReset').addEventListener('click', () => { session.reset(); drawInteractive(result); });
  }

  function siteTable(session) {
    const rows = session.siteReconstruction(10);
    const lines = ['Most-parsimonious reconstruction on the current tree', ''];
    lines.push('site   changes   states at the root');
    for (const r of rows) lines.push(String(r.site).padStart(4) + '   ' + String(r.changes).padStart(5) + '     ' + r.states.join(','));
    lines.push('');
    lines.push('Apply rearrangements above and watch this table change.');
    return lines.join('\n');
  }

  // ------------------------------------------------------------- downloads

  function download(name, content, type) {
    const url = URL.createObjectURL(new Blob([content], { type: type || 'text/plain' }));
    const a = document.createElement('a');
    a.href = url;
    a.download = name;
    a.click();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  }

  $('downloadNewick').addEventListener('click', () => {
    const tree = activeSession ? activeSession.tree : (last && last.kind === 'tree' ? last.tree : null);
    if (!tree) return;
    download(current.id.toLowerCase() + '.nwk', C.toNewick(tree) + '\n');
  });

  $('downloadText').addEventListener('click', () => {
    if (!last) return;
    download(current.id.toLowerCase() + '.txt', $('programText').textContent + '\n');
  });

  $('runProgram').addEventListener('click', run);
  $('programSearch').addEventListener('input', e => renderCatalogue(e.target.value));
  $('loadProgramExample').addEventListener('click', () => {
    if (!current) return;
    if (activeInput !== null) {
      if (store[activeInput] !== undefined && !EXAMPLES[current.input]) return;
      $('programInput').value = EXAMPLES[current.input] || EXAMPLES[activeInput];
    }
    if (current.needsTrees > 0) $('programTree').value = EXAMPLES[current.needsTrees > 1 ? 'trees' : 'tree'];
    $('programError').textContent = '';
  });

  ensureStatusNodes();
  renderCatalogue('');
  select('DNAPARS');
})();
