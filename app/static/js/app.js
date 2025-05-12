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
        // Get nested ReactionStep objects from dependencies
        const nestedSteps = Object.values(rx.dependencies)
            .filter(step => step !== null && step.reaction); // Ensure valid ReactionStep
        if (nestedSteps.length > 0) {
            const sub = buildTree(nestedSteps); // Recurse on dependencies
            li.appendChild(sub);
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

    results.forEach(item => {
        if (!item.chain || !item.product) {
            console.warn('Invalid result item:', item);
            return;
        }
        treeRoot.appendChild(buildTree(item.chain));
        const prodDiv = document.createElement('div');
        prodDiv.classList.add('product-card');
        prodDiv.innerHTML = `
            <h3>${item.product.chemical_name || 'Unknown'}</h3>
            <p>CAS: ${item.product.casid || 'Unknown'}</p>
            <p>Formula: ${item.product.formula || 'Unknown'}</p>
            <p>SMILES: ${item.product.smiles || 'Unknown'}</p>
        `;
        productList.appendChild(prodDiv);
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