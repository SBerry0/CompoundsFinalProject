from fastapi import FastAPI, Request
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from fastapi.responses import JSONResponse
from pydantic import BaseModel
from typing import List

from reactions import get_all_chemicals, find_paths

app = FastAPI()

app.mount("/static", StaticFiles(directory="static"), name="static")
templates = Jinja2Templates(directory="templates")


@app.get("/")
async def index(request: Request):
    return templates.TemplateResponse("index.html", {"request": request})


class Selection(BaseModel):
    chemicals: List[str]


@app.get("/api/chemicals")
async def list_chemicals():
    return get_all_chemicals()


@app.post("/api/paths")
async def compute_paths(sel: Selection):
    tree = find_paths(sel.chemicals)
    return JSONResponse(content=tree)


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)