import io

from fastapi import FastAPI, Request, HTTPException
from fastapi.responses import HTMLResponse, JSONResponse, PlainTextResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from typing import List
import json

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.AllChem import EmbedMolecule, UFFOptimizeMolecule

from app.models import ReactionChain
from app.chemicals import Chemical
# from app.reactions import Reaction
from app.reaction_logic import get_products

app = FastAPI()

# Mount static files folder
app.mount("/static", StaticFiles(directory="/Users/sohumberry/PycharmProjects/CompoundsFinalProject/app/static"), name="static")

# Templates directory
templates = Jinja2Templates(directory="app/templates")

# Load reactant data on startup
with open("/Users/sohumberry/PycharmProjects/CompoundsFinalProject/app/data/reactants.json") as f:
    reactant_list = [Chemical(**item) for item in json.load(f)]

@app.get("/", response_class=HTMLResponse)
async def index(request: Request):
    # We'll serve a mostly blank template; JS will populate reactants
    return templates.TemplateResponse("index.html", {"request": request})

@app.get("/reactants", response_class=JSONResponse)
async def get_reactants():
    # Return list of all possible reactants
    out = [r.to_dict() for r in reactant_list]
    # print("test: ", out)
    return JSONResponse([r.to_dict() for r in reactant_list])

@app.post("/run_reaction", response_class=JSONResponse)
async def run_reaction(reactants: List[str]):
    # reactants: list of CAS identifiers selected
    reactants_lower = [r.lower() for r in reactants]
    selected = [r for r in reactant_list if r.chemical_name.lower() in reactants_lower]

    print("vs......", reactants_lower)
    # Compute products and chains
    result: ReactionChain = get_products(selected)
    # Serialize output
    output = []
    for product, chain in result.items():
        output.append({
            "product": product.to_dict(),
            "chain": [rx.to_dict() for rx in chain]
        })
    print(output)
    return JSONResponse(output)

@app.get("/get_3d_model", response_class=PlainTextResponse)
async def get_3d_model(smiles: str):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string.")

        mol = Chem.AddHs(mol)
        if EmbedMolecule(mol) != 0:
            raise RuntimeError("3D coordinate embedding failed.")
        UFFOptimizeMolecule(mol)

        pdb_block = Chem.MolToPDBBlock(mol)
        return pdb_block
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"Failed to generate 3D model: {str(e)}")