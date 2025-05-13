import json
import pubchempy as pcp
import pandas as pd

from app.chemicals import Chemical


class Reaction:
    def __init__(self, id: str, reactants: [Chemical], products: [Chemical], reaction_type: str, conditions: str, reference: str):
        self.id = id
        self.reactants = reactants
        self.products = products
        self.reaction_type = reaction_type
        self.conditions = conditions
        self.reference = reference

    def __str__(self):
        str = ""
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

    def to_dict(self):
        return {'id': self.id,
                'reactants': [chem.to_dict() for chem in self.reactants],
                'products': [chem.to_dict() for chem in self.products],
                'reaction_type': self.reaction_type,
                'conditions': self.conditions,
                'reference': self.reference}

    @classmethod
    def from_dict(cls, d):
        return cls(d['id'], [Chemical.from_dict(c) for c in d['reactants']],
            [Chemical.from_dict(c) for c in d['products']], d['reaction_type'], d['conditions'], d['reference'])

def resolve_all(id, reactants_smiles, products_smiles, reaction_type, conditions, reference):
    # time.sleep(random.uniform(0.2, 1.0))
    return {
        'id': id,
        'reactants': [create_chemical(smile) for smile in reactants_smiles],
        'products': [create_chemical(smile) for smile in products_smiles],
        'reaction_type': reaction_type,
        'conditions': conditions,
        'reference': reference
    }

chemical_cache = {}

def create_chemical(smiles):
    if smiles in chemical_cache:
        return chemical_cache[smiles]

    try:
        compound = pcp.get_compounds(smiles, 'smiles')
        if compound:
            c = compound[0]
            chem = Chemical(
                c.cid,
                smiles,
                c.molecular_formula,
                c.synonyms[0] if c.synonyms else "No name"
            )
        else:
            chem = Chemical("ERROR", smiles, "ERROR", "ERROR")
    except Exception as e:
        print(f"Error resolving {smiles}: {e}")
        chem = Chemical("ERROR", smiles, "ERROR", "ERROR")
    # Dynamic Programming
    chemical_cache[smiles] = chem
    return chem

def clean_json(reaction):
    for key in ['reactants', 'products']:
        if key in reaction and isinstance(reaction[key], list):
            for compound in reaction[key]:
                if 'chemical_name' in compound and isinstance(compound['chemical_name'], str):
                    compound['chemical_name'] = compound['chemical_name'].title()
    return reaction

if __name__ == '__main__':
    # Read from input file
    with open('data/common_reactions.json', 'r', encoding='utf-8') as f:
        data = json.load(f)

    # Apply the transformation to each reaction (assuming list of reactions)
    if isinstance(data, list):
        updated_data = [clean_json(rxn) for rxn in data]
    else:
        updated_data = clean_json(data)

    # Write to output file
    with open('reactions_updated.json', 'w', encoding='utf-8') as f:
        json.dump(updated_data, f, indent=2, ensure_ascii=False)
    # reactions = []
    #
    # with open('raw_reactions.json', 'r') as file:
    #     data = json.load(file)
    #     print(len(data))
    #     for d in data:
    #         item = resolve_all(d['Reaction ID'], d['Reactants'], d['Products'], d['Reaction Type'], d['Conditions'], d['Reference'])
    #         reactions.append(Reaction(item['id'], item['reactants'], item['products'], item['reaction_type'], item['conditions'], item['reference']))
    #         print(item)
    #
    # with open("common_reactions.json", 'w') as f:
    #     json.dump([r.to_dict() for r in reactions], f, indent=2)



def pubchem():
    # Load your CSV
    compounds = pd.read_csv('disk/common_chemicals.csv')
    reactions = []
    reaction_count = 0
    for cas in compounds['CAS Registry Number']:
        try:
            compound = pcp.get_compounds(cas, 'cid')[0]
            reaction_data = compound.reactions  # Fetch reaction data
            for reaction in reaction_data:
                print(reaction)
                reactions.append({
                    'RXN_ID': f'RXN-{str(reaction_count).zfill(3)}',
                    'Reactants': reaction['reactants'],
                    'Products': reaction['products'],
                    'Type': reaction.get('type', 'Unknown'),
                    'Conditions': reaction.get('conditions', 'Not specified')
                })
                reaction_count += 1
        except:
            print(f"No data for CAS {cas}")

    # Save to CSV
    pd.DataFrame(reactions).to_csv('disk/full_reactions.csv', index=False)


def ord(compounds):
    # Load ORD JSON (replace with actual file path)
    with open('ord_data.json') as f:
        ord_data = json.load(f)

    reactions = []
    reaction_count = 101
    for reaction in ord_data['reactions']:
        if any(cas in reaction['reactants'] for cas in compounds['CAS Registry Number']):
            reactions.append({
                'Reaction ID': f'RXN-{str(reaction_count).zfill(3)}',
                'Reactants': reaction['reactants'],
                'Products': reaction['products'],
                'Reaction Type': reaction.get('type', 'Unknown'),
                'Conditions': reaction.get('conditions', 'Not specified'),
                'Reference': 'ORD'
            })
            reaction_count += 1

    pd.DataFrame(reactions).to_csv('ord_reactions.csv', index=False)