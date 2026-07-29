# Agent Instructions

This repository is the public paper-companion snapshot for
`topological-css-code-decoupling`.

## Boundaries

- Keep the repository history clean; do not import private Git history.
- Do not copy paper sources into this repository.
- Do not add private URLs, local absolute paths, serialized caches, secrets,
  debug dumps, or unresolved publication placeholders.
- Keep generated reproduction output under ignored build directories.

## Commands

Use focused checks while editing:

```bash
python3 -m unittest discover -s tests/release -v
python3 tools/verify_release.py
```

Python/SageMath checks run from `python/` after installing `.[test]` in Sage.
Julia/Oscar checks run with `julia --project=julia/ToricBuilder`.

## Escalation

If a verification command cannot be run because SageMath, Julia, or a locked
dependency is unavailable, report that limitation with the exact command and
do not claim the check passed.
