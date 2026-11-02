const example = `>Human
ATGCTAGCTAGCTACGATCGATCGTACGATCG
>Chimpanzee
ATGCTAGCTAGCTACGATCGATCGTACGATCA
>Gorilla
ATGCTAGCTAGCTACGATCGATCATACGATCA
>Orangutan
ATGCTAGCTAGCTACGATCAATCATACGATCA
>Gibbon
ATGCTAGCTAGTTACGATCAATCATACGATCA`;
const $ = id => document.getElementById(id);
let result = null;
$('sequences').value = example;
function updateStatus() {
  const count = ($('sequences').value.match(/^>/gm) || []).length;
  $('inputStatus').textContent = count ? `${count} sequence${count === 1 ? '' : 's'} in editor` : 'Ready for your sequences';
  // Results reflect the last run, not subsequent edits.
  if (result) { result = null; $('download').disabled = true; $('treeArea').innerHTML = '<div class="empty"><div class="empty-icon">✳</div><strong>Ready for a fresh tree</strong><p>Your input changed. Build the tree again<br>to see updated results.</p></div>'; ['statSeq','statSites','statMethod'].forEach(id => $(id).textContent = '—'); }
}
$('sequences').addEventListener('input', updateStatus);
$('model').addEventListener('change', updateStatus);
$('loadExample').addEventListener('click', () => { $('sequences').value = example; $('error').textContent = ''; updateStatus(); });
$('file').addEventListener('change', async event => {
  const file = event.target.files[0]; if (!file) return;
  try { $('sequences').value = await file.text(); $('error').textContent = ''; updateStatus(); } catch (_) { $('error').textContent = 'Could not read that file.'; }
  event.target.value = '';
});
function drawTree(tree, count) {
  const ns = 'http://www.w3.org/2000/svg';
  const svg = document.createElementNS(ns, 'svg');
  const width = 650, height = Math.max(300, count * 35 + 40), left = 35, right = 430;
  svg.setAttribute('viewBox', `0 0 ${width} ${height}`); svg.setAttribute('role', 'img'); svg.setAttribute('aria-label', 'UPGMA phylogenetic tree');
  svg.style.minHeight = height + 'px';
  let leaf = 0;
  const maxHeight = tree.height;
  function xFor(node, depth) { return maxHeight > 0 ? left + (right - left) * (1 - node.height / maxHeight) : left + Math.min(depth * 45, right - left); }
  function line(x1,y1,x2,y2) { const el = document.createElementNS(ns,'line'); for (const [k,v] of Object.entries({x1,y1,x2,y2,stroke:'#86a98a','stroke-width':1.7})) el.setAttribute(k,v); svg.appendChild(el); }
  function visit(node, depth) {
    const x = xFor(node, depth);
    if (!node.children) {
      const y = 25 + leaf++ * (height - 50) / Math.max(count - 1, 1);
      const dot = document.createElementNS(ns,'circle'); dot.setAttribute('cx',x); dot.setAttribute('cy',y); dot.setAttribute('r',4); dot.setAttribute('fill','#699b70'); svg.appendChild(dot);
      const label = document.createElementNS(ns,'text'); label.setAttribute('x',x+10); label.setAttribute('y',y+4); label.setAttribute('fill','#456c50'); label.setAttribute('font-size','12'); label.setAttribute('font-family','DM Sans, sans-serif'); label.textContent = node.name.length > 27 ? node.name.slice(0, 26) + '…' : node.name; const title = document.createElementNS(ns,'title'); title.textContent = node.name; label.appendChild(title); svg.appendChild(label); return {x,y};
    }
    const children = node.children.map(c => visit(c, depth + 1)); const y = (children[0].y + children[1].y) / 2;
    line(x, children[0].y, x, children[1].y); children.forEach(c => line(x,c.y,c.x,c.y));
    const dot = document.createElementNS(ns,'circle'); dot.setAttribute('cx',x); dot.setAttribute('cy',y); dot.setAttribute('r',3); dot.setAttribute('fill','#c0d2bd'); svg.appendChild(dot); return {x,y};
  }
  visit(tree, 0); $('treeArea').replaceChildren(svg);
}
$('run').addEventListener('click', () => {
  try {
    const next = Phylo.analyze($('sequences').value, $('model').value);
    result = next; drawTree(result.tree, result.records.length);
    $('statSeq').textContent = result.records.length; $('statSites').textContent = result.sites; $('statMethod').textContent = 'UPGMA';
    $('download').disabled = false; $('error').textContent = '';
  } catch (error) { updateStatus(); $('error').textContent = error.message; }
});
$('download').addEventListener('click', () => {
  if (!result) return;
  const url = URL.createObjectURL(new Blob([Phylo.newick(result.tree) + '\n'], {type:'text/plain'}));
  const a = document.createElement('a'); a.href = url; a.download = 'phoenix-tree.nwk'; a.click(); setTimeout(() => URL.revokeObjectURL(url), 1000);
});
updateStatus();
