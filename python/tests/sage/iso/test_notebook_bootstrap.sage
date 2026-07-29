import _bootstrap  # noqa: F401

import json
from pathlib import Path

repo_root = Path(_bootstrap.REPO_ROOT)
notebook_root = repo_root / "results"

missing_bootstrap = []
for notebook_path in sorted(notebook_root.rglob("*.ipynb")):
    if ".ipynb_checkpoints" in notebook_path.parts:
        continue
    notebook = json.loads(notebook_path.read_text())
    code_cells = [
        "".join(cell.get("source", []))
        for cell in notebook.get("cells", [])
        if cell.get("cell_type") == "code"
    ]
    if not any("isomorphism" in source for source in code_cells):
        continue
    first_code = code_cells[0] if code_cells else ""
    if "Sage-Decoder notebook import bootstrap" not in first_code:
        missing_bootstrap.append(str(notebook_path.relative_to(repo_root)))

assert missing_bootstrap == [], missing_bootstrap

print("test_notebook_bootstrap: ok")
