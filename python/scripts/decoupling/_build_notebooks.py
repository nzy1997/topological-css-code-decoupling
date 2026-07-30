"""Build the three explanatory notebooks from shared examples and results."""

from pathlib import Path

import nbformat


SOURCE_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_DIRECTORY = SOURCE_ROOT / "results" / "decoupling"


def markdown(text):
    """Return one Markdown notebook cell."""
    return nbformat.v4.new_markdown_cell(text)


def code(text):
    """Return one Sage/Python notebook cell."""
    return nbformat.v4.new_code_cell(text)


def notebook(cells):
    """Return a SageMath notebook with deterministic metadata."""
    document = nbformat.v4.new_notebook(cells=cells)
    document.metadata.kernelspec = {
        "display_name": "SageMath 10.7",
        "language": "sage",
        "name": "sagemath-10.7",
    }
    document.metadata.language_info = {
        "codemirror_mode": {"name": "ipython", "version": 3},
        "file_extension": ".py",
        "mimetype": "text/x-python",
        "name": "python",
        "nbconvert_exporter": "python",
        "pygments_lexer": "ipython3",
        "version": "3.13.3",
    }
    return document


def color_code_notebook(label):
    """Build one color-code notebook that consumes its canonical result JSON."""
    constant = f"COLOR_CODE_{label}"
    title = "6.6.6" if label == "666" else "4.8.8"
    return notebook(
        [
            markdown(
                f"# {title} color-code decoupling\n\n"
                "This notebook reads the shared paper example and the canonical "
                "reproduction output. Run the corresponding script before executing "
                "the notebook."
            ),
            code(
                "from pathlib import Path\n"
                "import json\n"
                "import sys\n"
                "from sage.all import Matrix\n"
                "SOURCE_ROOT = next(\n"
                "    candidate\n"
                "    for candidate in (Path.cwd(), *Path.cwd().parents)\n"
                "    if (candidate / 'decoupling' / '__init__.py').is_file()\n"
                ")\n"
                "sys.path.insert(0, str(SOURCE_ROOT))\n"
                "for module_name in tuple(sys.modules):\n"
                "    if module_name == 'decoupling' or "
                "module_name.startswith('decoupling.'):\n"
                "        del sys.modules[module_name]\n"
                f"from decoupling.paper_examples import {constant}\n"
                "from decoupling import R, x, y\n"
                "from sage.misc.sage_eval import sage_eval\n\n"
                "def matrix_from_strings(rows):\n"
                "    return Matrix(\n"
                "        R,\n"
                "        [[R(sage_eval(entry, locals={'x': x, 'y': y})) "
                "for entry in row] for row in rows],\n"
                "    )\n\n"
                f"case = {constant}\n"
                f"result_path = SOURCE_ROOT / 'results' / 'decoupling' / "
                f"'color_code_{label}.json'\n"
                "payload = json.loads(result_path.read_text(encoding='utf-8'))\n"
                "record = payload['result']\n"
                "assert record['verification_status'] == 'passed'\n"
                "assert record['superlattice_basis'] == "
                "[list(vector) for vector in case.superlattice_basis]\n"
                "record['actual_q'], record['actual_p_x'], "
                "record['actual_p_z'], record['actual_t']"
            ),
            markdown("## Paper input and paper-specified superlattice"),
            code(
                "input_excitation_map = case.excitation_map()\n"
                "input_excitation_map, case.superlattice_basis, "
                "record['actual_square_period']"
            ),
            markdown("## Computed standard complex"),
            code(
                "h_x_tilde = matrix_from_strings(record['h_x_tilde'])\n"
                "h_z_tilde_dagger = matrix_from_strings("
                "record['h_z_tilde_dagger'])\n"
                "h_x_tilde, h_z_tilde_dagger"
            ),
            markdown("## Computed inverse chain isomorphisms"),
            code(
                "psi_2_inverse = matrix_from_strings(record['psi_2_inverse'])\n"
                "psi_1_inverse = matrix_from_strings(record['psi_1_inverse'])\n"
                "psi_0_inverse = matrix_from_strings(record['psi_0_inverse'])\n"
                "psi_2_inverse, psi_1_inverse, psi_0_inverse"
            ),
            markdown(
                "The chain isomorphisms are non-unique. Different valid Gröbner "
                "bases, quotient bases, and pivot choices can produce matrices that "
                "do not match the paper entry by entry. The reproduction script "
                "therefore verifies the forward and inverse chain equations, mutual "
                "invertibility, and symplectic/QCA decoupling identities."
            ),
            code(
                "record['actual_psi_1_inverse_degree'], "
                "record['decoupling_runtime_seconds'], "
                "record['verification_status']"
            ),
        ]
    )


def bb_notebook():
    """Build the notebook that consumes the complete 59-row BB result."""
    return notebook(
        [
            markdown(
                "# BB-code decoupling results\n\n"
                "The complete 59-row inventory lives in "
                "`decoupling.paper_examples`; this notebook does not duplicate any "
                "Laurent polynomial or table value. The canonical result applies "
                "the documented benchmark-19 runtime cutoff."
            ),
            code(
                "from pathlib import Path\n"
                "import json\n"
                "import sys\n"
                "SOURCE_ROOT = next(\n"
                "    candidate\n"
                "    for candidate in (Path.cwd(), *Path.cwd().parents)\n"
                "    if (candidate / 'decoupling' / '__init__.py').is_file()\n"
                ")\n"
                "sys.path.insert(0, str(SOURCE_ROOT))\n"
                "for module_name in tuple(sys.modules):\n"
                "    if module_name == 'decoupling' or "
                "module_name.startswith('decoupling.'):\n"
                "        del sys.modules[module_name]\n"
                "from decoupling.paper_examples import BB_CODE_INSTANCES\n\n"
                "result_path = SOURCE_ROOT / 'results' / 'decoupling' / 'bb_codes.json'\n"
                "payload = json.loads(result_path.read_text(encoding='utf-8'))\n"
                "records = payload['rows']\n"
                "assert len(BB_CODE_INSTANCES) == 59\n"
                "assert payload['metadata']['catalog_row_count'] == 59\n"
                "assert len(records) == payload['metadata']['row_count']\n"
                "assert payload['metadata']['runtime_cutoff_reference_id'] "
                "== 'benchmark-19'\n"
                "assert payload['metadata']['complete']\n"
                "assert payload['metadata']['all_passed']\n"
                "len(BB_CODE_INSTANCES), len(records), "
                "payload['metadata']['excluded_longer_runtime_row_ids']"
            ),
            markdown("## Shared paper inputs"),
            code(
                "[case.as_reference_record() for case in BB_CODE_INSTANCES[:3]]"
            ),
            markdown("## Actual implementation outputs"),
            code(
                "columns = (\n"
                "    'id', 'actual_square_period', 'superlattice_basis', "
                "'actual_q',\n"
                "    'actual_p_x', 'actual_p_z', 'actual_t',\n"
                "    'actual_psi_1_inverse_degree', "
                "'decoupling_runtime_seconds',\n"
                "    'verification_status',\n"
                ")\n"
                "[{column: row[column] for column in columns} for row in records[:10]]"
            ),
            markdown(
                "Degree and runtime are measured outputs, not equality oracles. "
                "Chain isomorphisms are non-unique, so the validator checks the "
                "paper's structural values and mathematical identities instead of "
                "requiring entry-by-entry agreement."
            ),
            code(
                "assert all(row['verification_status'] == 'passed' for row in records)\n"
                "max(row['actual_q'] for row in records), "
                "max(row['actual_psi_1_inverse_degree'] for row in records)"
            ),
        ]
    )


def write_notebook(path, document):
    """Write one notebook with stable indentation and UTF-8 text."""
    path.write_text(nbformat.writes(document, version=4), encoding="utf-8")


write_notebook(OUTPUT_DIRECTORY / "666_color_code.ipynb", color_code_notebook("666"))
write_notebook(OUTPUT_DIRECTORY / "488_color_code.ipynb", color_code_notebook("488"))
write_notebook(OUTPUT_DIRECTORY / "bb_code_instances.ipynb", bb_notebook())
