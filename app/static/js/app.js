async function fetchReactants() {
    const res = await fetch('/reactants');
    if (!res.ok) {
        console.error('Failed to fetch reactants:', res.status);
        return [];
    }
    return res.json();
}

function createButton(chem) {
    const btn = document.createElement('button');
    btn.textContent = chem.chemical_name;
    btn.classList.add('reactant-btn');

    if (chem.chemical_name) {
        btn.dataset.chemicalName = chem.chemical_name;
    } else {
        console.warn("Missing chemical name for chemical:", chem);
    }

    btn.addEventListener('click', () => btn.classList.toggle('selected'));
    return btn;
}

function renderReactants(list) {
    const container = document.getElementById('reactant-container');
    if (!container) {
        console.error('Reactant container not found');
        return;
    }
    list.forEach(chem => {
        container.appendChild(createButton(chem));
    });
}

function buildTree(chain) {
    const ul = document.createElement('ul');
    chain.forEach(rx => {
        const li = document.createElement('li');
        li.textContent = `${rx.reaction.id}: ${rx.reaction.reaction_type} (${rx.reaction.conditions})`;

        // Create a nested <ul> for reaction details and dependencies
        const detailsUl = document.createElement('ul');

        // Build pretty reaction string: reactants + " → " + products
        const reactionLi = document.createElement('li');
        const reactantStr = rx.reaction.reactants && Array.isArray(rx.reaction.reactants)
            ? rx.reaction.reactants.map(r => r.formula || 'Unknown').join(' + ')
            : 'None';
        const productStr = rx.reaction.products && Array.isArray(rx.reaction.products)
            ? rx.reaction.products.map(p => p.formula || 'Unknown').join(' + ')
            : 'None';
        reactionLi.textContent = `${reactantStr} → ${productStr}`;
        detailsUl.appendChild(reactionLi);

        // Add Reactants label and names
        if (rx.reaction.reactants && Array.isArray(rx.reaction.reactants) && rx.reaction.reactants.length > 0) {
            const reactantsLabelLi = document.createElement('li');
            reactantsLabelLi.textContent = 'Reactants:';
            reactantsLabelLi.style.fontWeight = 'bold';
            detailsUl.appendChild(reactantsLabelLi);

            rx.reaction.reactants.forEach(r => {
                const nameLi = document.createElement('li');
                nameLi.textContent = r.chemical_name || 'Unknown';
                detailsUl.appendChild(nameLi);
            });
        }

        // Add Products label and names
        if (rx.reaction.products && Array.isArray(rx.reaction.products) && rx.reaction.products.length > 0) {
            const productsLabelLi = document.createElement('li');
            productsLabelLi.textContent = 'Products:';
            productsLabelLi.style.fontWeight = 'bold';
            detailsUl.appendChild(productsLabelLi);

            rx.reaction.products.forEach(p => {
                const nameLi = document.createElement('li');
                nameLi.textContent = p.chemical_name || 'Unknown';
                detailsUl.appendChild(nameLi);
            });
        }

        // Get nested ReactionStep objects from dependencies
        const nestedSteps = Object.values(rx.dependencies)
            .filter(step => step !== null && step.reaction);
        if (nestedSteps.length > 0) {
            const sub = buildTree(nestedSteps);
            detailsUl.appendChild(sub);
        }

        // Append the details <ul> if it has children
        if (detailsUl.children.length > 0) {
            li.appendChild(detailsUl);
        }

        ul.appendChild(li);
    });
    return ul;
}

function displayResults(results) {
    const treeRoot = document.getElementById('reaction-tree');
    const productList = document.getElementById('product-list');
    if (!treeRoot || !productList) {
        console.error('Reaction tree or product list container not found');
        return;
    }

    treeRoot.innerHTML = '';
    productList.innerHTML = '';

    if (!Array.isArray(results)) {
        console.error('Invalid results format:', results);
        return;
    }

    results.forEach((item, index) => {
        if (!item.chain || !item.product) {
            console.warn('Invalid result item:', item);
            return;
        }
        treeRoot.appendChild(buildTree(item.chain));
        const prodDiv = document.createElement('div');
        prodDiv.classList.add('product-card');
        // Add unique ID for 3D viewer
        const viewerId = `viewer-${index}`;
        prodDiv.innerHTML = `
            <h3>${item.product.chemical_name || 'Unknown'}</h3>
            <p>CAS: ${item.product.casid || 'Unknown'}</p>
            <p>Formula: ${item.product.formula || 'Unknown'}</p>
            <p>SMILES: ${item.product.smiles || 'Unknown'}</p>
            <div id="${viewerId}" class="viewer-3d" style="width: 200px; height: 200px;"></div>
        `;
        productList.appendChild(prodDiv);

        // Load 3D model if SMILES is available
        if (item.product.smiles) {
            fetch(`/get_3d_model?smiles=${encodeURIComponent(item.product.smiles)}`)
                .then(res => {
                    if (!res.ok) {
                        throw new Error(`Failed to fetch 3D model: ${res.status}`);
                    }
                    return res.text(); // Get PDB data as text
                })
                .then(pdbData => {
                    const viewer = $3Dmol.createViewer(viewerId);
                    viewer.addModel(pdbData, 'pdb');
                    viewer.setStyle({}, { stick: {} });
                    viewer.zoomTo();
                    viewer.render();
                })
                .catch(err => {
                    console.error('Error loading 3D model:', err);
                    prodDiv.querySelector(`#${viewerId}`).innerHTML = 'Failed to load 3D model';
                });
        } else {
            prodDiv.querySelector(`#${viewerId}`).innerHTML = 'No SMILES data';
        }
    });
}

async function runReaction() {
    const selected = Array.from(document.querySelectorAll('.reactant-btn.selected'))
        .map(btn => btn.dataset.chemicalName)
        .filter(name => name);
    if (selected.length === 0) {
        console.warn('No reactants selected');
        return;
    }
    console.log('Selected reactants:', JSON.stringify(selected));
    try {
        const res = await fetch('/run_reaction', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(selected)
        });
        if (!res.ok) {
            console.error('Failed to run reaction:', res.status, await res.text());
            return;
        }
        const data = await res.json();
        console.log('Response data:', data);
        displayResults(data);
    } catch (error) {
        console.error('Error running reaction:', error);
    }
}

document.getElementById('run-btn')?.addEventListener('click', runReaction);

// Initialize
fetchReactants().then(renderReactants).catch(error => {
    console.error('Error initializing reactants:', error);
});