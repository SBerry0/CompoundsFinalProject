import itertools
import json
from typing import List, Set
from reactions import Reaction
from chemicals import Chemical

class ReactionStep:
    def __init__(self, reaction: Reaction, dependencies: dict):
        self.reaction = reaction
        self.dependencies = dependencies  # dict: reactant → ReactionStep or None

def find_products(available_reactants: Set[str], reactions: List[Reaction]):
    known_chemicals = set(available_reactants)
    chemical_history = {} # chemical: (reaction, reactants)
    print("Finding products....")

    new_found = True
    while new_found:
        new_found = False
        for rxn in reactions:
            if rxn.reactants.issubset(known_chemicals):
                for product in rxn.products:
                    if product not in known_chemicals:
                        known_chemicals.add(product)
                        chemical_history[product] = (rxn, rxn.reactants.copy())
                        new_found = True
    print("Products found.")

    def get_chain(product: str):
        # If the given product is in the available reactants, it's an original material and should be excluded
        if product in available_reactants:
            return None
        reaction, reactants = chemical_history[product]
        # Recursively find everything that was required to create this product
        dependencies = {r: get_chain(r) for r in reactants}
        return ReactionStep(reaction.id, reactants, product, dependencies)

    return { product: get_chain(product) for product in chemical_history }


def print_reaction_chain(step: ReactionStep, indent: int = 0):
    if step is None:
        return
    pad = '  ' * indent
    print(f"{pad}{step.product} <= via reaction {step.reaction_id} using {step.reactants}")
    for reactant, substep in step.dependencies.items():
        print_reaction_chain(substep, indent + 1)


def get_all_chemicals():
    with open('data/common_chemicals.json', 'r') as f:
        data = json.load(f)
    return [Chemical.from_dict(d) for d in data]


if __name__ == "__main__":
    chemicals = get_all_chemicals()

    with open('data/common_reactions.json', 'r') as f:
        data = json.load(f)
    reactions = [Reaction.from_dict(d) for d in data]

    starting = set(itertools.islice(chemicals, 5))

    chains = find_products(starting, reactions)
    print(chains)
    for product, chain in chains.items():
        print(f"\nChain for {product}:")
        print_reaction_chain(chain)
