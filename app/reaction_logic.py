import itertools
import json
from typing import List, Set
from app.chemicals import Chemical
from app.reactions import Reaction

class ReactionStep:
    def __init__(self, reaction: Reaction, dependencies: dict):
        self.reaction = reaction
        self.dependencies = dependencies  # dict: reactant → ReactionStep or None

    def __iter__(self):
        # Yield the current ReactionStep
        yield self
        # Recursively yield from dependencies that are ReactionStep instances
        for step in self.dependencies.values():
            if isinstance(step, ReactionStep):
                yield from step

    def to_dict(self):
        # Convert dependencies, using chemical_name as keys
        dependencies_dict = {}
        for chemical, step in self.dependencies.items():
            key = chemical.chemical_name if isinstance(chemical, Chemical) else str(chemical)
            dependencies_dict[key] = step.to_dict() if isinstance(step, ReactionStep) else None

        reaction_dict = (
            self.reaction.to_dict()
            if hasattr(self.reaction, "to_dict")
            else {"id": getattr(self.reaction, "id", None)}
        )

        return {
            "reaction": reaction_dict,
            "dependencies": dependencies_dict
        }

def find_products(available_reactants: Set[Chemical], reactions: List[Reaction]):
    known_chemicals = set(available_reactants)
    chemical_history = {} # chemical: (reaction, reactants)
    print("Finding products....")

    new_found = True
    while new_found:
        new_found = False
        for rxn in reactions:
            if set(rxn.reactants).issubset(known_chemicals):
                for product in rxn.products:
                    if product not in known_chemicals:
                        known_chemicals.add(product)
                        chemical_history[product] = (rxn, rxn.reactants.copy())
                        new_found = True
    print("Products found.")

    def get_chain(product: Chemical):
        # If the given product is in the available reactants, it's an original material and should be excluded
        if product in available_reactants:
            return None
        reaction, reactants = chemical_history[product]
        # Recursively find everything that was required to create this product
        dependencies = {r: get_chain(r) for r in reactants}
        return ReactionStep(reaction, dependencies)

    return { product: get_chain(product) for product in chemical_history }


def print_reaction_chain(step: ReactionStep, indent: int = 0):
    if step is None:
        return
    pad = '  ' * indent
    print(f"{pad}{step.reaction.products} <= via {step.reaction.reaction_type} reaction {step.reaction.id} using {step.reaction.reactants}")
    for reactant, substep in step.dependencies.items():
        print_reaction_chain(substep, indent + 1)

def get_products(reactants: [Chemical]):
    with open('/Users/sohumberry/PycharmProjects/CompoundsFinalProject/app/data/common_reactions.json', 'r') as f:
        data = json.load(f)
    reactions = [Reaction.from_dict(d) for d in data]
    return find_products(set(reactants), reactions)

def get_all_chemicals():
    with open('data/reactants.json', 'r') as f:
        data = json.load(f)
    return [Chemical.from_dict(d) for d in data]


if __name__ == "__main__":
    chemicals = get_all_chemicals()

    with open('data/common_reactions.json', 'r') as f:
        data = json.load(f)
    reactions = [Reaction.from_dict(d) for d in data]

    # starting = set(itertools.islice(chemicals, 5))

    starting = set([c for c in chemicals if c.chemical_name in ["Acetic acid", "Ethanol"]])
    chains = get_products(list(starting))
    # chains = get_products(chemicals)
    print(chains)
    for product, chain in chains.items():
        print(f"\nChain for {product}:")
        print_reaction_chain(chain)
