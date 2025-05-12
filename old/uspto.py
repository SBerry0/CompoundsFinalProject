from itertools import product
from rxnmapper import RXNMapper
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw, rdChemReactions
import os

pd.set_option('display.max_rows', 100)  # Show up to 10 rows
pd.set_option('display.max_columns', 15)  # Show up to 5 columns
pd.set_option('display.width', 1000)  # Set display width to 100 characters
pd.set_option('display.max_colwidth', 500) # Set max column width to 20 characters


def full_rxn_data():
    set_a = pd.read_csv("../data/dataSetA.csv")
    set_b = pd.read_csv("../data/dataSetB.csv")
    set_1000 = pd.read_csv("../data/uspto_1k_TPL_train_valid.tsv", delimiter="\t")

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
    df_uspto = pd.read_csv("../data/dataSetB.csv")
    # print(df_ocr)
    # reactants_unclean = df_ocr.iloc[:,0]
    # products_unclean = df_ocr.iloc[:,1]
    print(df_uspto)

    rxns = df_uspto["rxnSmiles_Mapping_NameRxn"]
    is_completes = df_uspto["NameRxn_Mapping_Complete"]
    rxn_ids = df_uspto["patentID"]
    reaction_smiles = []
    ids = []
    for i in range(len(rxns)):
        if is_completes[i]:
            reaction_smiles.append(rxns[i])
            ids.append(rxn_ids[i])

    products = []
    reactants = []

    for reaction in reaction_smiles:
        rxn = rdChemReactions.ReactionFromSmarts(reaction, useSmiles=True)
        rxn_reactants = [Chem.MolToSmiles(rxn.GetReactantTemplate(i)) for i in range(rxn.GetNumReactantTemplates())]
        rxn_products = [Chem.MolToSmiles(rxn.GetProductTemplate(i)) for i in range(rxn.GetNumProductTemplates())]
        reactants.append(rxn_reactants)
        products.append(rxn_products)


    results = pd.DataFrame({
        "patent_id": ids,
        "reactants": reactants,
        "products": products,
        "reaction_smiles": reaction_smiles
    })
    results.to_csv("data/uspto_data.csv", index=False)


def inspect_data(data_path):
    data = pd.read_csv(data_path)

    print(data.shape)
    print(data)

# def analyze_reactions(SMILES_string):
#     reaction = rdChemReactions.ReactionFromSmarts(SMILES_string, useSmiles=True)
#     print(type(reaction))
#
#     # rxn = rdChemReactions.ReactionFromSmarts("[C:1](=O)[Cl:2].[NH2:3]>>[C:1](=O)[NH:3]", useSmiles=True)
#     mol1 = Chem.MolFromSmiles("CCCCCC")
#     mol2 = Chem.MolFromSmiles("CO")
#     mol3 = Chem.MolFromSmiles("C")
#     mol4 = Chem.MolFromSmiles("[ClH:28]")
#     mol5 = Chem.MolFromSmiles("[Pd]")
#     mol6 = Chem.MolFromSmiles("[ClH:28]")
#
#     # products = reaction.RunReactants((mol1, mol2, mol3, mol4, mol5, mol6))
#     products = reaction.RunReactant(Chem.MolFromSmiles("C"), 0)
#     print(products)
#     for product_set in products:
#         for product in product_set:
#             print("Product:", Chem.MolToSmiles(product))


# def map_reactions(data):
#     mapper = RXNMapper()
#     unmapped = pd.read_csv(data)
#
#     results = mapper.get_attention_guided_atom_maps(rxns=unmapped["reaction_smiles"], detailed_output=True)
#     # results = mapper.get_attention_guided_atom_maps(rxns=["CC(C)(C)NN.O=C1CCCCC1.[C-]#N>>CC(C)(C)NNC1(C#N)CCCCC1"])
#     mapped = []
#     # mapped_rxn = results[0]['mapped_rxn']
#     # print("Mapped Reaction:", mapped_rxn)
#     print(len(results))
#     for i in range(len(results)):
#         if i % 1000 == 0:
#             print(i, "/", len(results))
#         mapped.append(results[i]['mapped_rxn'])
#
#     output = pd.DataFrame(mapped)
#     output.to_csv("data/mapped.csv")


def count_heavy_atoms(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return mol.GetNumHeavyAtoms() if mol else 999  # 999 ensures error ones are excluded

def is_reaction_small_enough(reaction_smiles, max_heavy_atoms):
    try:
        parts = reaction_smiles.strip().split(">>")
        reactants = parts[0].split(".") if parts[0] else []
        products = parts[2].split(".") if len(parts) > 2 else []

        all_smiles = reactants
        return all(count_heavy_atoms(sm) <= max_heavy_atoms for sm in all_smiles)
    except Exception as e:
        return False

def filter_reactions(input_file, output_file, mol_size_threshold):
    input = pd.read_csv(input_file)
    clean_smiles = []
    clean_reactants = []
    clean_products = []
    clean_ids = []
    for i in range(len(input["reaction_smiles"])):
        if is_reaction_small_enough(input["reaction_smiles"][i], mol_size_threshold):
            clean_smiles.append(input["reaction_smiles"][i])
            clean_reactants.append(input["reactants"][i])
            clean_products.append(input["products"][i])
            clean_ids.append((input["patent_id"][i]))

    results = pd.DataFrame({
        "patent_id": clean_ids,
        "reactants": clean_reactants,
        "products": clean_products,
        "reaction_smiles": clean_smiles
    })
    results.to_csv(output_file, index=False)

    print(f"Filtered reactions saved to: {output_file}")
    print(input.shape, "vs", results.shape)



if __name__ == "__main__":
    filter_reactions("../data/uspto_data.csv", "data/uspto_clean_data_20.csv", 20)
    # inspect_data("data/uspto_data.csv")
    # inspect_data("data/uspto_clean_data.csv")
    # analyze_reactions("C.CCCCCC.CO.O=C(OCc1ccccc1)[NH:1][CH2:2][CH2:3][CH2:4][CH2:5][C@@H:6]([C:7]([O:8][CH3:9])=[O:10])[NH:11][C:12](=[O:13])[NH:14][c:15]1[cH:16][c:17]([O:18][CH3:19])[cH:20][c:21]([C:22]([CH3:23])([CH3:24])[CH3:25])[c:26]1[OH:27].[ClH:28].[Pd]>>[ClH:28].[NH2:1][CH2:2][CH2:3][CH2:4][CH2:5][C@@H:6]([C:7]([O:8][CH3:9])=[O:10])[NH:11][C:12](=[O:13])[NH:14][c:15]1[cH:16][c:17]([O:18][CH3:19])[cH:20][c:21]([C:22]([CH3:23])([CH3:24])[CH3:25])[c:26]1[OH:27]")
    # map_reactions("data/full_ocr_data")
    # inspect_data("data/full_ocr_data")
    # split_data()
    # analyze_rxn()
    # full_rxn_data()