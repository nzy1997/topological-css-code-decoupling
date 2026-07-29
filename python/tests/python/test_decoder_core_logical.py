"""Check logical classification helpers that do not require Sage."""

from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[2]
repo_root = str(REPO_ROOT)
if repo_root not in sys.path:
    sys.path.insert(0, repo_root)

from decoder_core.logical import classify_verified_attempt  # noqa: E402


class FakeMatrix:
    def __mul__(self, correction):
        return correction.syndrome


class FakeClassifier:
    def is_logical_failure(self, residual):
        return bool(residual)


class FakeError:
    def __add__(self, correction):
        return correction.residual


class FakeCorrection:
    def __init__(self, syndrome, residual):
        self.syndrome = syndrome
        self.residual = residual


logical_failure, decode_failure = classify_verified_attempt(
    FakeMatrix(),
    FakeClassifier(),
    FakeError(),
    "expected",
    FakeCorrection("expected", True),
    decoder_name="unitary_decouple",
)
assert logical_failure
assert not decode_failure

try:
    classify_verified_attempt(
        FakeMatrix(),
        FakeClassifier(),
        FakeError(),
        "expected",
        FakeCorrection("wrong", False),
        decoder_name="unitary_decouple",
    )
except ValueError as exc:
    message = str(exc)
    assert "unitary_decouple" in message
    assert "syndrome" in message
else:
    raise AssertionError("Expected a verified decoder syndrome mismatch to raise.")

print("test_decoder_core_logical: ok")
