import ast
import json

def load_smiles_from_json(json_file):
    with open(json_file, 'r') as f:
        data = json.load(f)

    out = set(item["smiles"].lower() for item in data if "smiles" in item)
    out.add("co")
    return out


import csv

import re

# Strip atom mapping from SMILES
def strip_atom_mapping(smiles):
    # Remove things like ":1", ":23", etc. from atoms
    smiles = smiles[1:-1]
    smiles = smiles.replace("'", "")
    smiles = re.sub(r'\:([0-9]+)', '', smiles)
    return smiles.strip().lower()


def filter_reactions(csv_file, smiles_set, output_file, reactant_col='reactants', product_col='products', delimiter=', ', partial=True):
    total_rows = 0
    kept_rows = 0

    with open(csv_file, 'r', newline='') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        writer = csv.DictWriter(outfile, fieldnames=reader.fieldnames)
        writer.writeheader()

        for row in reader:
            total_rows += 1
            reactants = {strip_atom_mapping(s) for s in row[reactant_col].split(delimiter) if s.strip()}
            # reactants = set(ast.literal_eval(row['reactants']))
            # if total_rows == 1:
                # print(reactants, smiles_set)
            if partial:
                if reactants & smiles_set:
                    writer.writerow(row)
                    kept_rows += 1
            else:
                if reactants.issubset(smiles_set):
                    writer.writerow(row)
                    kept_rows += 1


    print(f"Total reactions in input: {total_rows}")
    print(f"Reactions matching partial SMILES: {kept_rows}")



smiles_set = load_smiles_from_json('data/common_chemicals.json')
filter_reactions('../data/filteredRxnSize/uspto_clean_data_30.csv', smiles_set, 'filtered_reactions_r_full.csv', partial=False)
