"""Reproduce the paper's 6.6.6 color-code decoupling example."""

from pathlib import Path
import sys


SOURCE_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_DIRECTORY = Path(__file__).resolve().parent
for path in (SOURCE_ROOT, SCRIPT_DIRECTORY):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from decoupling.paper_examples import COLOR_CODE_666  # noqa: E402
from _color_code_runner import run_color_code  # noqa: E402


record = run_color_code(COLOR_CODE_666, SOURCE_ROOT)
print(
    f"6.6.6 color-code reproduction: {record['verification_status']} "
    f"in {record['actual_runtime_seconds']:.6f} s"
)
