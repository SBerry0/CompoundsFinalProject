import json
import csv

import cirpy
import pubchempy as pcp

class Chemical:
    def __init__(self, casid: str, smiles: str, formula: str, chemical_name: str):
        self.casid = casid
        self.smiles = smiles
        self.formula = formula
        self.chemical_name = chemical_name

    def __eq__(self, other):
        if isinstance(other, Chemical):
            return self.chemical_name.lower() == other.chemical_name.lower()
        return NotImplemented

    def __hash__(self):
        return hash(self.chemical_name.lower())

    def __ne__(self, other):
        return not self.__eq__(other)

    def __repr__(self):
        return f"{self.chemical_name} ({self.casid} - {self.smiles})"

    def __str__(self):
        return self.smiles

    def to_dict(self):
        return {'casid': self.casid,
                'smiles': self.smiles,
                'formula': self.formula,
                'chemical_name': self.chemical_name}

    @classmethod
    def from_dict(cls, d):
        return cls(d['casid'], d['smiles'], d['formula'], d['chemical_name'])


def resolve_all(casid, name, formula):
    # time.sleep(random.uniform(0.2, 1.0))
    return {
        'casid': casid,
        'smiles': cirpy.resolve(casid, 'smiles'),
        'formula': formula,
        'chemical_name': name,
    }


if __name__ == '__main__':
    chemicals = []


    with open('disk/common_chemicals.csv', 'r') as file:
        reader = csv.reader(file)
        next(reader)
        for row in reader:
            item = resolve_all(row[2], row[0], row[1])
            chemicals.append(Chemical(item['casid'], item['smiles'], item['formula'], item['chemical_name']))

            print(row)

    with open("data/common_chemicals.json", 'w') as f:
        json.dump([c.to_dict() for c in chemicals], f, indent=2)

# with open('chemicals.json', 'r') as f:
#     data = json.load(f)
#
# loaded_compounds = [Chemical.from_dict(d) for d in data]