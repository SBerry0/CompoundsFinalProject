import os

import pandas as pd
import pubchempy as pcp
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
import time

pd.set_option('display.max_columns', None)


list_file = "../data/cpdat_v4.0/cpdat_v4.0_list_presence_data.csv"
functional_file = "../data/cpdat_v4.0/cpdat_v4.0_functional_use_data.csv"
output_file = "../data/cpdat_merged_files.csv"

checkpoint_interval = 2000
checkpoint_file_prefix = "chemicals_checkpoint_"
final_output = "chemicals_large.json"


# Load both files
list_presence_df = pd.read_csv(list_file, low_memory=False)
print(list_presence_df.shape)
list_presence_df = list_presence_df.dropna(subset=['curated_casrn'])



# Build the CID → name map
ids = {row['curated_casrn']: row['curated_chemical_name'] for _, row in list_presence_df.iterrows()}


class Chemical:
    def __init__(self, cas: str, smiles: str, formula: str, chemical_name: str, names: [str]):
        self.casid = cas
        self.smiles = smiles
        self.formula = formula
        self.chemical_name = chemical_name
        self.names = names

    def to_dict(self):
        return {
            'casid': self.casid,
            'smiles': self.smiles,
            'formula': self.formula,
            'chemical_name': self.chemical_name,
            'name': self.names
        }

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

# Retry logic included for robustness
def resolve_pubchem(cid, retries=3, delay=4):
    for attempt in range(retries):
        try:
            compound = pcp.Compound.from_cid(cid)
            return {
                'cid': cid,
                'smiles': compound.isomeric_smiles,
                'name': compound.iupac_name,
                'formula': compound.molecular_formula
            }
        except Exception as e:
            if attempt < retries - 1:
                time.sleep(delay)
            else:
                return {
                    'cid': cid,
                    'smiles': None,
                    'name': None,
                    'formula': None,
                    'error': str(e)
                }

print(resolve_pubchem(next(iter(ids))))
# Run in parallel
# chemicals = []
# with ThreadPoolExecutor(max_workers=20) as executor:
#     futures = {executor.submit(resolve_pubchem, cid): cid for cid in ids.keys()}
#
#     for i, future in enumerate(as_completed(futures)):
#         if i % 200 == 0:
#             print(f"{i}/{len(futures)} complete")
#         cid = futures[future]
#         try:
#             res = future.result()
#             chemicals.append(Chemical(
#                 cid, res['smiles'], res['formula'], ids[cid], res['name']
#             ))
#         except Exception as e:
#             chemicals.append(Chemical(cid, None, None, ids[cid], [f"Error: {e}"]))

# Save to JSON
# with open('chemicals_large.json', 'w') as f:
#     json.dump([c.to_dict() for c in chemicals], f, indent=2)
