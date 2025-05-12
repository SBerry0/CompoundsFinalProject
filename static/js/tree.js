async function loadChemicals() {
  const resp = await fetch('/api/chemicals');
  return resp.json();
}

function createSelector(options) {
  const sel = document.createElement('select');
  options.forEach(ch => {
    const o = document.createElement('option'); o.value = o.text = ch;
    sel.append(o);
  });
  return sel;
}

(async () => {
  const chems = await loadChemicals();
  const container = document.getElementById('selectors');
  document.getElementById('add-btn').onclick = () => container.append(createSelector(chems));
  // Add two selectors by default
  container.append(createSelector(chems));
  container.append(createSelector(chems));

  document.getElementById('go-btn').onclick = async () => {
    const selected = Array.from(container.querySelectorAll('select')).map(s => s.value);
    const res = await fetch('/api/paths', {
      method: 'POST', headers: {'Content-Type':'application/json'},
      body: JSON.stringify({ chemicals: selected })
    });
    const treeData = await res.json();
    renderTree(treeData);
    listProducts(treeData);
  };
})();

function renderTree(data) {
  document.getElementById('tree').innerHTML = '';
  const root = d3.hierarchy({ children: data }, d => d.children);
  const treeLayout = d3.tree().size([400, 200]);
  treeLayout(root);

  const svg = d3.select('#tree')
    .append('svg')
      .attr('width', 500)
      .attr('height', 300)
    .append('g')
      .attr('transform', 'translate(40,0)');

  // Links
  svg.selectAll('line')
    .data(root.links())
    .enter()
    .append('line')
      .attr('x1', d => d.source.x)
      .attr('y1', d => d.source.y)
      .attr('x2', d => d.target.x)
      .attr('y2', d => d.target.y);

  // Nodes
  svg.selectAll('circle')
    .data(root.descendants())
    .enter()
    .append('circle')
      .attr('cx', d => d.x)
      .attr('cy', d => d.y)
      .attr('r', 5);

  // Labels with reaction names
  svg.selectAll('text')
    .data(root.descendants().slice(1))
    .enter()
    .append('text')
      .attr('x', d => (d.parent.x + d.x) / 2)
      .attr('y', d => (d.parent.y + d.y) / 2)
      .text(d => d.data.reaction);
}

function listProducts(tree) {
  const list = document.getElementById('product-list');
  list.innerHTML = '';
  function traverse(nodes, depth = 1) {
    nodes.forEach(n => {
      const li = document.createElement('li');
      li.textContent = `${n.name} (steps: ${depth})`;
      list.append(li);
      if (n.children) traverse(n.children, depth+1);
    });
  }
  traverse(tree);
}