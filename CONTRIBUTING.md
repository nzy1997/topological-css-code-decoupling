# Contributing

This repository is a paper companion, so changes should keep mathematical
claims, generated artifacts, and documentation traceable to explicit commands.

## Setup

Install the relevant environment before making changes:

```bash
cd python
sage -pip install -e '.[test]'
```

```bash
julia --project=julia/ToricBuilder -e 'using Pkg; Pkg.instantiate()'
```

## Validation

Run release checks before proposing changes:

```bash
python3 -m unittest discover -s tests/release -v
python3 tools/verify_release.py
```

Run the Python/SageMath and Julia/Oscar tests that cover your change. If a
result file changes, document the command that produced it and keep generated
bulk output under ignored build directories.

## Scope

Do not add private development notes, machine-local paths, serialized caches,
debug dumps, credentials, or paper-source files. Prefer small pull requests
with a clear reproduction or verification command.
