import json
import os

import pandas as pd
pd.set_option('display.max_columns', None)


list_file = "../data/cpdat_v4.0/cpdat_v4.0_list_presence_data.csv"
functional_file = "../data/cpdat_v4.0/cpdat_v4.0_functional_use_data.csv"
output_file = "../data/cpdat_merged_files.csv"

checkpoint_interval = 2000
checkpoint_file_prefix = "chemicals_checkpoint_"
final_output = "chemicals_large2.json"


# Load both files
# list_presence_df = pd.read_csv(list_file, low_memory=False)
# print(list_presence_df.shape)
# list_presence_df = list_presence_df.dropna(subset=['curated_casrn'])
# print(list_presence_df.shape)
functional_df = pd.read_csv(functional_file, low_memory=False)
functional_df = functional_df.dropna(subset=['curated_casrn'])

# Select relevant columns
functional_cols = functional_df[[
    'data_source', 'data_group', 'curated_chemical_name', 'curated_casrn', 'dtxsid',
    'raw_functional_use', 'function_category', 'component_name'
]]

# list_cols = list_presence_df[[
#     'data_source', 'curated_chemical_name', 'curated_casrn', 'dtxsid',
#     'data_group', 'organization', 'component_name', 'keyword_set'
# ]]

# for col in list_cols:
#     print(list_presence_df.iloc[9][col])

# list_presence_df = list_presence_df.dropna(subset=['component_name'])
# print(list_presence_df.shape)



# for col in list_cols:
#     print(list_presence_df.iloc[9][col])

# print(list_presence_df.head(5))
# list_clean = list_presence_df[['curated_chemical_name', 'dtxsid', 'data_']]


# for col in list_presence_df.columns:
#     unique_vals = list_presence_df[col].unique()
#     print(f"{col}: {unique_vals}")



from collections import defaultdict


ids = {}

for _, row in functional_df.iterrows():
    ids[row['curated_casrn']] = row['function_category']


# print(ids)
print(len(ids))

import pubchempy as pcp

def resolve_pubchem(key):
    try:
        compound = pcp.Compound.from_cid(key)
        return {
            'cid': key,
            'smiles': compound.isomeric_smiles,
            'name': compound.iupac_name,
            'formula': compound.molecular_formula
        }
    except IndexError:
        return {'cid': key, 'smiles': None, 'name': None, 'formula': None}
    except Exception as e:
        return {'cid': key, 'smiles': None, 'name': None, 'formula': None, 'error': str(e)}



import cirpy
from concurrent.futures import ThreadPoolExecutor, as_completed
import time
import random


def resolve_all(key):
    # time.sleep(random.uniform(0.2, 1.0))
    return {
        'cid': key,
        'smiles': cirpy.resolve(key, 'smiles'),
        'name': cirpy.resolve(key, 'names'),
        'formula': cirpy.resolve(key, 'formula')
    }

# results = []
# i=0
# with ThreadPoolExecutor(max_workers=10) as executor:
#     for result in executor.map(resolve_all, ids.keys()):
#         if i % 200 == 0:
#             print(i, "/", len(ids))
#         results.append(result)
#         i+=1


# for key, value in ids.items():
#     if i % 200 == 0:
#         print(i, "/", len(ids))
#     smiles = cirpy.resolve(key, 'smiles')
#     name = cirpy.resolve(key, 'names')
#     formula = cirpy.resolve(key, 'formula')
#     # print(name)
#     ids[key].append(smiles)
#     ids[key].append(name)
#     ids[key].append(formula)
#     i+=1

# print(ids)


class Chemical:
    def __init__(self, cas: str, smiles: str, formula: str, chemical_name: str, names: [str]):
        self.casid = cas
        self.smiles = smiles
        self.formula = formula
        self.chemical_name = chemical_name
        self.names = names

    def to_dict(self):
        return {'casid': self.casid,
                'smiles': self.smiles,
                'formula': self.formula,
                'chemical_name': self.chemical_name,
                'name': self.names}

    @classmethod
    def from_dict(cls, d):
        return cls(d['casid'], d['smiles'], d['formula'], d['chemical_name'], d['name'])


def load_last_checkpoint():
    files = sorted([
        f for f in os.listdir() if f.startswith(checkpoint_file_prefix) and f.endswith(".json")
    ], key=lambda x: int(x.split("_")[-1].split(".")[0]))
    if not files:
        return [], set()
    last_file = files[-1]
    with open(last_file, 'r') as f:
        data = json.load(f)
        chemicals = [Chemical.from_dict(d) for d in data]
        processed_ids = {chem.casid for chem in chemicals}
        print(f"Resuming from {last_file} with {len(chemicals)} entries")
        return chemicals, processed_ids


chemicals, processed = load_last_checkpoint()

# 🧪 Prepare remaining IDs
remaining = [cid for cid in ids.keys() if cid not in processed]

# 🚀 Run in parallel
with ThreadPoolExecutor(max_workers=10) as executor:
    futures = {executor.submit(resolve_all, cid): cid for cid in remaining}

    for i, future in enumerate(as_completed(futures), start=len(chemicals)):
        cid = futures[future]
        try:
            res = future.result()
            chemical = Chemical(cid, res['smiles'], res['formula'], ids[cid], res['name'])
        except Exception as e:
            chemical = Chemical(cid, None, None, ids[cid], [f"Error: {e}"])

        chemicals.append(chemical)

        if (i + 1) % checkpoint_interval == 0:
            ckpt_file = f"{checkpoint_file_prefix}{i + 1}.json"
            print(f"Saving checkpoint to {ckpt_file}...")
            with open(ckpt_file, 'w') as f:
                json.dump([c.to_dict() for c in chemicals], f, indent=2)

# 🧾 Final save
with open(final_output, 'w') as f:
    json.dump([c.to_dict() for c in chemicals], f, indent=2)

print(f"✅ Finished! Total chemicals saved: {len(chemicals)}")


# chemicals = []
# # results = [resolve_pubchem(key) for key in ids.keys()]
# results = []
# for i, key in enumerate(ids.keys()):
#
#     # results.append(resolve_pubchem(key))
#     res = resolve_all(key)
#     results.append(res)
#     if i % 200 == 0:
#         print(i, "/", len(ids))
#         print(res)
# for res in results:
#     chemicals.append(Chemical(res['cid'], res['smiles'], res['formula'], ids[res['key']], res['name']))
#
# import json
#
# with open('chemicals_large.json', 'w') as f:
#     json.dump([c.to_dict() for c in chemicals], f, indent=2)

# with open('chemicals.json', 'r') as f:
#     data = json.load(f)
#
# loaded_compounds = [Chemical.from_dict(d) for d in data]