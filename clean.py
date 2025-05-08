import pandas as pd
pd.set_option('display.max_columns', None)


list_file = "data/cpdat_v4.0/cpdat_v4.0_list_presence_data.csv"
functional_file = "data/cpdat_v4.0/cpdat_v4.0_functional_use_data.csv"
output_file = "data/cpdat_merged_files.csv"
# Load both files
list_presence_df = pd.read_csv(list_file, low_memory=False)
print(list_presence_df.shape)
list_presence_df = list_presence_df.dropna(subset=['dtxsid'])
print(list_presence_df.shape)
# functional_df = pd.read_csv(functional_file, low_memory=False)

# Select relevant columns
# functional_cols = functional_df[[
#     'data_source', 'data_group', 'curated_chemical_name', 'curated_casrn', 'dtxsid',
#     'raw_functional_use', 'function_category', 'component_name'
# ]]

list_cols = list_presence_df[[
    'data_source', 'curated_chemical_name', 'curated_casrn', 'dtxsid',
    'data_group', 'organization', 'component_name', 'keyword_set'
]]

for col in list_cols:
    print(list_presence_df.iloc[9][col])

list_presence_df = list_presence_df.dropna(subset=['component_name'])
print(list_presence_df.shape)



for col in list_cols:
    print(list_presence_df.iloc[9][col])

print(list_presence_df.head(5))

# list_clean = list_presence_df[['curated_chemical_name', 'dtxsid', 'data_']]

print('--------------------------------------------------------')
print('\n\n')


for col in list_presence_df.columns:
    unique_vals = list_presence_df[col].unique()
    print(f"{col}: {unique_vals}")



from collections import defaultdict


ids = defaultdict(list)

for _, row in list_presence_df.iterrows():
    ids[row['curated_casrn']] = [row['curated_chemical_name'], row["component_name"]]


print(ids)
print(len(ids))


import cirpy

for key, value in ids.items():
    smiles = cirpy.resolve(key, 'smiles')
    name = cirpy.resolve(key, 'names')
    # print(name)
    ids[key].append(smiles)
    ids[key].append(name)

print(ids)


class Chemical:
    def __init__(self, cas: str, smiles: str, chemical_name: str, component_name: str, names: [str]):
        self.casid = cas
        self.smiles = smiles
        self.chemical_name = chemical_name
        self.component_name = component_name
        self.names = names

    def to_dict(self):
        return {'casid': self.casid,
                'smiles': self.smiles,
                'chemical_name': self.chemical_name,
                'component_name': self.component_name,
                'name': self.names}

    @classmethod
    def from_dict(cls, d):
        return cls(d['casid'], d['smiles'], d['chemical_name'], d['component_name'], d['name'])

chemicals = []

for key, value in ids.items():
    chemicals.append(Chemical(key, value[2], value[0], value[1], value[3]))

import json

with open('chemicals.json', 'w') as f:
    json.dump([c.to_dict() for c in chemicals], f, indent=2)

# with open('chemicals.json', 'r') as f:
#     data = json.load(f)
#
# loaded_compounds = [Chemical.from_dict(d) for d in data]