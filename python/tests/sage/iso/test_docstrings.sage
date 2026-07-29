import _bootstrap  # noqa: F401

import ast
from pathlib import Path


repo_root = Path(__file__).resolve().parents[3]
missing_docstrings = []

for path in sorted((repo_root / "isomorphism").rglob("*.py")):
    tree = ast.parse(path.read_text())
    for node in ast.walk(tree):
        if isinstance(node, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            if ast.get_docstring(node) is None:
                missing_docstrings.append(
                    f"{path.relative_to(repo_root)}:{node.lineno}:{node.name}"
                )

assert not missing_docstrings, "Missing docstrings: " + ", ".join(missing_docstrings)

print("test_docstrings: ok")
