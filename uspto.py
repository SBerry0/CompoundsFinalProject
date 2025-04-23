from itertools import product
from rxnmapper import RXNMapper
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw, rdChemReactions


pd.set_option('display.max_rows', 100)  # Show up to 10 rows
pd.set_option('display.max_columns', 15)  # Show up to 5 columns
pd.set_option('display.width', 1000)  # Set display width to 100 characters
pd.set_option('display.max_colwidth', 500) # Set max column width to 20 characters


def full_rxn_data():
    set_a = pd.read_csv("data/dataSetA.csv")
    set_b = pd.read_csv("data/dataSetB.csv")
    set_1000 = pd.read_csv("data/uspto_1k_TPL_train_valid.tsv", delimiter="\t")

    # Don't seem super useful
    # iim_b = pd.read_csv("data/incompleteIndigoMappings_dataSetB.csv")
    # iim_b_results = pd.read_csv("data/incompleteIndigoMappings_dataSetB_results.csv")

    # print(set_a)
    # print("____________________________________________________________________")
    # print(set_b)

    # size of 683 rxns
    rxn_smiles_a = set_a["rxn_Smiles"]
    # size of 50,000 rxns
    rxn_smiles_b = set_b["rxnSmiles_Mapping_NameRxn"]

    rxn_smiles_1k = set_1000["original_rxn"]

    print(type(rxn_smiles_a))
    print(rxn_smiles_a)
    print(rxn_smiles_b)
    print(rxn_smiles_1k)

    sample_smiles = rxn_smiles_1k.iloc[10]

    print("sample: ", sample_smiles)
    print(type(sample_smiles))
    # rxn = Chem.MolFromSmiles(sample_smiles)

    rxn = rdChemReactions.ReactionFromSmarts(sample_smiles, useSmiles=True)
    return analyze_rxn(rxn)


def analyze_rxn(rxn):
    # Print the number of reactants/products
    print(f"Reactants: {rxn.GetNumReactantTemplates()}")
    print(f"Products: {rxn.GetNumProductTemplates()}")
    reactants = []
    # Print SMILES of individual molecules
    print("\nReactants:")
    for i in range(rxn.GetNumReactantTemplates()):
        reactants.append(Chem.MolToSmiles(rxn.GetReactantTemplate(i)))
        print(Chem.MolToSmiles(rxn.GetReactantTemplate(i)))
    products = []
    print("\nProducts:")
    for i in range(rxn.GetNumProductTemplates()):
        products.append(Chem.MolToSmiles(rxn.GetProductTemplate(i)))
        print(Chem.MolToSmiles(rxn.GetProductTemplate(i)))
    img = Draw.ReactionToImage(rxn)
    img.show()
    return reactants, products


def split_data():
    df_ocr = pd.read_csv("data/ocrtrain.csv")
    # print(df_ocr)
    reactants_unclean = df_ocr.iloc[:,0]
    products_unclean = df_ocr.iloc[:,1]

    reactants = [s.replace(" ", "") for s in reactants_unclean]
    products = [s.replace(" ", "") for s in products_unclean]

    reaction_smiles = [""] * len(reactants)
    reactions = [rdChemReactions.ChemicalReaction] * len(reactants)
    for i in range(len(reactants)):
        if i % 5000 == 0:
            print(i, "/", len(reactants))
        reaction_smiles[i] = reactants[i] + ">>" + products[i]
        # print(reaction_smiles[i])
        reactions[i] = rdChemReactions.ReactionFromSmarts(reaction_smiles[i], useSmiles=True)

    results = pd.DataFrame({
        "reactants": reactants,
        "products": products,
        "reaction_smiles": reaction_smiles
    })
    results.to_csv("data/full_ocr_data", index=False)

    for i in range(20):
        analyze_rxn(reactions[i])


def inspect_data(data_path):
    data = pd.read_csv(data_path)

    print(data.shape)
    print(data)

def analyze_reactions(SMILES_string):
    reaction = rdChemReactions.ReactionFromSmarts(SMILES_string, useSmiles=True)
    print(type(reaction))

    # rxn = rdChemReactions.ReactionFromSmarts("[C:1](=O)[Cl:2].[NH2:3]>>[C:1](=O)[NH:3]", useSmiles=True)
    mol1 = Chem.MolFromSmiles("OCCCl")
    mol2 = Chem.MolFromSmiles("CC(C)CS(=O)(=O)Cl")
    products = reaction.RunReactants((mol1, mol2))

    for product_set in products:
        for product in product_set:
            print("Product:", Chem.MolToSmiles(product))


def map_reactions(data):
    mapper = RXNMapper()
    unmapped = pd.read_csv(data)

    results = mapper.get_attention_guided_atom_maps(rxns=unmapped["reaction_smiles"], detailed_output=True)
    # results = mapper.get_attention_guided_atom_maps(rxns=["CC(C)(C)NN.O=C1CCCCC1.[C-]#N>>CC(C)(C)NNC1(C#N)CCCCC1"])
    mapped = []
    # mapped_rxn = results[0]['mapped_rxn']
    # print("Mapped Reaction:", mapped_rxn)
    print(len(results))
    for i in range(len(results)):
        if i % 1000 == 0:
            print(i, "/", len(results))
        mapped.append(results[i]['mapped_rxn'])

    output = pd.DataFrame(mapped)
    output.to_csv("data/mapped.csv")

if __name__ == "__main__":
    # analyze_reactions("CC(C)CS(=O)(=O)Cl.OCCCl>>CC(C)CS(=O)(=O)OCCCl")
    map_reactions("data/full_ocr_data")
    # inspect_data("data/full_ocr_data")
    # split_data()
    # analyze_rxn()
    # full_rxn_data()