"""Check tracked isomorphism notebooks use the public supplement API."""

import _bootstrap  # noqa: F401

import json
from pathlib import Path


repo_root = Path(_bootstrap.REPO_ROOT)
notebook_root = repo_root / "results" / "isomorphism"
notebook_paths = sorted(notebook_root.glob("*.ipynb"))

assert [path.name for path in notebook_paths] == [
    "488_color_code.ipynb",
    "666_color_code.ipynb",
    "bb_code_instances.ipynb",
]

legacy_fragments = (
    ".total_" + "maps",
    ".q" + "_inverse",
    '"q' + '_shape"',
    '"q' + '_inverse_computed"',
    ".input_coarse_matrix",
    ".clearing",
    ".post_inverse_maps",
    ".inverse_maps.hx_standard",
    ".inverse_maps.hz_standard",
    ".maps.hx_standard",
    ".maps.hz_standard",
)

for notebook_path in notebook_paths:
    notebook = json.loads(notebook_path.read_text())
    code_cells = [
        "".join(cell.get("source", []))
        for cell in notebook.get("cells", [])
        if cell.get("cell_type") == "code"
    ]
    for cell_index, code in enumerate(code_cells):
        compile(code, f"{notebook_path.name}:code-cell-{cell_index}", "exec")
    source = "\n".join(
        "".join(cell.get("source", []))
        for cell in notebook.get("cells", [])
    )
    for fragment in legacy_fragments:
        assert fragment not in source, (notebook_path.name, fragment)
    assert ".inverse_maps" in source, notebook_path.name
    assert "phi1_inverse" in source, notebook_path.name

print("test_notebook_api: ok")
