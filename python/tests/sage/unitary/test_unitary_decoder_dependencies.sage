"""Ensure unitary decoder is independent of matching decoder packages."""

import _bootstrap


repo_root = _bootstrap.REPO_ROOT
package_root = repo_root / "unitary_decouple_decoder"
legacy_prefix = "bp" + "_matching" + "_decoder"

for path in package_root.glob("*.py"):
    source = path.read_text(encoding="utf-8")
    for forbidden in (
        legacy_prefix + "_old_version",
        legacy_prefix,
        "matching" + "_decoder",
    ):
        assert forbidden not in source

print("test_unitary_decoder_dependencies: ok")
