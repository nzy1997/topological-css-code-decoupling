# Python Upstream Synchronization Design

## Goal

Bring the public `topological-css-code-decoupling` repository's Python/Sage
implementation to the final state of private-repository commit `cc23bae`
(`Align Python decoupling implementation with paper`). The public repository
will expose only the new API; the superseded package names and reproduction
entry points will not be retained as compatibility aliases.

## Source and target

- Source: `topological-code-decoupling/python_code/` at commit `cc23bae` (the
  current source tree includes this commit through merge commit `8436df0`).
- Target: `topological-css-code-decoupling/python/`.
- The source Python tree is authoritative for implementation code,
  reproduction scripts, notebooks, and tracked Python result artifacts.
- The public repository remains authoritative for repository-level license,
  confirmed authorship, citation metadata, CI, contribution guidance, and
  release verification infrastructure.

## Synchronization rules

1. Replace the `isomorphism` package with `decoupling`.
2. Replace `unitary_decouple_decoder` with
   `unitary_decouple_based_decoder`.
3. Synchronize `decoder_core`, the new paper-facing APIs, reproduction
   scripts, notebooks, and tracked result layout from the source tree.
4. Apply source-side deletions: remove superseded Python tests, scripts,
   notebooks, and result files instead of retaining stale copies.
5. Adapt package metadata for the public repository: keep the confirmed
   software author and MIT release status while adopting the new package
   names, version, description, dependencies, and package discovery rules.
6. Do not copy private-tree release staging files into `python/`, including
   `.decoupling-public-release`, `CITATION.cff`, `MANIFEST.sha256`, and
   `RELEASE_BLOCKERS.md`.
7. Update public-repository integration files whose commands or paths become
   stale, including the root README, reproduction documentation, Makefile,
   CI workflow, release tests, and verification tooling where required.
8. Regenerate the public root `MANIFEST.sha256` after the tree reaches its
   final state.

## Public API outcome

The published Python packages will be:

- `decoupling`
- `unitary_decouple_based_decoder`
- `decoder_core`

Imports from `isomorphism` and `unitary_decouple_decoder` are intentionally
removed. The API terminology follows the paper: forward maps are `psi`,
inverse maps are `psi_inverse`, and the main result is
`DecouplingUnitaryResult`.

## Verification

The synchronized repository must pass:

```bash
python3 -m unittest discover -s tests/release -v
python3 tools/verify_release.py
```

Run focused pure-Python checks and source-provided Sage reproduction or smoke
checks when their runtimes are available. If SageMath or optional dependencies
are unavailable, record the exact unrun command and do not report it as
passing. Before completion, inspect the final diff for private URLs, absolute
local paths, caches, secrets, debug output, publication placeholders, and
unintentional changes outside the Python integration surface.

## Commit structure

Keep the design approval separate from implementation. The implementation
commit will contain the synchronized Python tree, required public integration
updates, regenerated manifest, and no unrelated Julia or paper-source changes.
