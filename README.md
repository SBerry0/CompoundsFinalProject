# Reaction Explorer

Reaction Explorer is my Advanced Topics in Computer Science final project. Select however many reactants from the selection field, then click **Run Reaction**.

The program will display all reactions that will take place (in the dataset) and through which process, even including multi-step reactions where products of one reaction is used as a reactant for another.

Finally, all the created products will be displayed with a 3D rendering of the molecular structure if available.

## Usage

This program uses FastAPI to run on localhost. As such, the program is run with the command `uvicorn main:app`. Then navigate to http://127.0.0.1:8000/ to view the program.
(Ignore the old folder, that's for Mr. Blick)

## Libraries

- PubChemPy
- RDKit
- CIRPy
- FastAPI
- Pandas

## Data

Final:
- Web Scraping

Attempts:
- USPTO
- OCR
- BioBricks
- CPDat
