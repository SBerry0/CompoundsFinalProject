from collections import defaultdict
from typing import List, Set
import ast
import pandas as pd


class Reaction:
    def __init__(self, id: str, reactants: List[str], products: List[str]):
        # Turn the list into a set to remove duplicates and take advantage of Python methods
        self.id = id
        self.reactants: Set[str] = set(reactants)
        self.products: Set[str] = set(products)

    #this is so dumb
    def __str__(self):
        str = "Reaction id: " + self.id + "\n"
        for i in range(len(self.reactants)):
            str += list(self.reactants)[i]
            if i+1 < len(self.reactants):
                str += " + "
        str += " >> "
        for i in range(len(self.products)):
            str += list(self.products)[i]
            if i+1 < len(self.products):
                str += " + "
        return str

class ReactionStep:
    def __init__(self, reaction_id: str, reactants: Set[str], product: str, dependencies: dict):
        self.reaction_id = reaction_id
        self.reactants = reactants
        self.product = product
        self.dependencies = dependencies  # dict: reactant → ReactionStep or None

# find the products that can also be created with products from other reactions (no keeping track of reactions)
def find_multi_step_products(available_reactants: Set[str], available_reactions: List[Reaction]):
        known_chemicals = set(available_reactants)
        new_found = True

        while new_found:
            new_found = False
            for rxn in available_reactions:
                if rxn.reactants.issubset(known_chemicals):
                    added = rxn.products - known_chemicals
                    if added:
                        known_chemicals.update(added)
                        new_found = True

        return known_chemicals - available_reactants  # Only return new products


def find_products(available_reactants: Set[str], available_reactions: List[Reaction]):
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


# for testing from chatGPT
def print_reaction_chain(step: ReactionStep, indent: int = 0):
    if step is None:
        return
    pad = '  ' * indent
    print(f"{pad}{step.product} <= via reaction {step.reaction_id} using {step.reactants}")
    for reactant, substep in step.dependencies.items():
        print_reaction_chain(substep, indent + 1)

sample_reactions = [
    # Reaction("US08978554B2", ['Cl[C:1]([CH:2]=[CH2:3])=[O:4]', '[NH2:5][CH2:6][C:7](=[O:8])[OH:9]', '[Na+]', '[OH-]'], ['[C:1]([CH:2]=[CH2:3])(=[O:4])[NH:5][CH2:6][C:7](=[O:8])[OH:9]'])
    Reaction("1", ["C(C1C(C(C(C(O1)O)O)O)O)O"], ["CCO.CCO", "O=C=O.O=C=O"]),
    Reaction("2", ["Cl", "N"], ["[NH4+].[Cl-]"]),
    Reaction("3", ["NaHCO3", "CC(=O)O"], ["CC(=O)[O-].[Na+]", "O", "O=C=O"]),
    Reaction("4", ["Cu", "N", "O", "O=O"], ["[Cu(NH3)4](OH)2"]),
    Reaction("5", ["NaHCO3"], ["Na2CO3", "O", "O=C=O"]),
    Reaction("6", ["O=C=O", "[H][H]"], ["C", "O"])
]

def make_reactions(filepath):
    print("Setting up reactions....")
    print("Using file \"" + filepath + "\"....")
    if not filepath:
        return sample_reactions

    df = pd.read_csv(filepath)
    reactions = []
    for row in df.itertuples():
        # ast.literal_eval from ChatGPT
        reactions.append(Reaction(row.patent_id, ast.literal_eval(row.reactants), ast.literal_eval(row.products)))
    print(reactions[1])
    print("Reactions are set.\n")
    return reactions


if __name__ == "__main__":
    starting = {"Cl", "N", "NaHCO3", "[H][H]"}
    # for rxn in reactions:
    #     for reactant in rxn.reactants:
    #         print(reactant, end=", ")

    # print(find_multi_step_products(starting, reactions))

    # paths = find_products({"Cl", "N", "NaHCO3"}, reactions)
    # for product, path_info in paths.items():
    #     print(f"{product} was created via reaction {path_info['reaction_id']} using {path_info['reactants']}")

    reactions = make_reactions("data/filteredReactantSize/uspto_clean_data_10.csv")
    print("Running reactions with starting chemicals", starting)
    chains = find_products(starting, reactions)
    print(chains)
    for product, chain in chains.items():
        print(f"\nChain for {product}:")
        print_reaction_chain(chain)

