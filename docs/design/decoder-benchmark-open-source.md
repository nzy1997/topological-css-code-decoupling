# Decoder benchmark open-source design

Status: implemented as archival releases `paper-decoder-2026.1` and
`paper-decoder-benchmark-2026.2`. This document preserves the design
rationale; the authoritative commands and current limitations are in
`docs/reproduction.md`.

## Need

The article's decoder figure is currently reproducible only at the plotting
layer. The companion repository contains the plotted arrays, but it does not
publish a runnable version of the historical Julia decoder and benchmark
pipeline that produced them. This leaves the benchmark protocol, confidence
intervals, matching weights, timing scope, software revisions, and raw-data
selection insufficiently documented.

The release must let an external researcher do three distinct things without
access to private repositories or local paths:

1. identify the exact historical implementation and inputs used for the
   article;
2. run a small deterministic end-to-end smoke benchmark; and
3. run the full sample-to-failures protocol and compare the new statistical
   result with the archived article data.

The historical run did not record random seeds. The release therefore does
not promise bit-for-bit regeneration of the archived Monte Carlo counts. It
does promise a fully specified, newly seeded protocol and preservation of the
original raw results with checksums and provenance.

## Repository responsibilities

The release uses three repositories with non-overlapping responsibilities.

### TensorQEC.jl

TensorQEC owns the decoder implementation used by the article:

- the historical Julia unitary-decouple decoder;
- the packed binary matrix operations used by the timed decoder path;
- the unweighted sparse-blossom matching adapter;
- decoder compilation and correctness smoke tests.

Create `paper/topological-css-decoder` from commit `9ffb8a9`. This point
contains the optimized decoder introduced by `18c43e6` and the test inputs
added immediately before the article's logical-error data were generated.
The later `BundleMatchingDecoder` development on `zy/basic_to_toric` is not
part of the historical implementation and is excluded from this release.

After the branch is reproducible on a clean machine, tag it
`paper-decoder-2026.1`. The companion and benchmark repositories must depend
on the immutable tag's commit, not on the branch name.

### DecoderBenchmarks.jl

DecoderBenchmarks owns experiment orchestration and statistical processing:

- construction of the article's BB-code family at toric-code distances 4, 6,
  8, and 10;
- sample-to-failures Monte Carlo execution;
- BP-OSD and unitary-decouple benchmark entrypoints;
- decoder-kernel timing;
- profile-likelihood confidence intervals;
- raw result selection and plot-data export.

Create `paper/topological-css-decoder` from commit `7a7d3bd`, the last commit
before the archived logical-error data were generated. Do not base the branch
on the later `zy/polish` worktree because that line mixes the historical
unitary-decouple experiment with later bundle-decoder experiments.

After the branch is reproducible and its selected raw data are verified, tag
it `paper-decoder-benchmark-2026.2`.

### topological-css-code-decoupling

The companion repository owns the paper-facing interface:

- exact external commit identifiers;
- the normalized article dataset;
- the provenance manifest and checksums;
- the reproduction guide and Make targets;
- regression checks that the normalized dataset still corresponds to the
  selected archived raw files;
- the final figure renderer.

It does not copy the TensorQEC or DecoderBenchmarks source trees. This avoids
duplicate implementations and preserves a clear authorship and maintenance
boundary.

## Historical benchmark contract

### Code family

The benchmark covers the BB code with

```text
H_X = (1 + x + x^-1*y, 1 + y + x*y)
```

at toric-code distances 4, 6, 8, and 10. Each exported curve must record its
physical qubit count and `[[n,k,d]]` parameters in addition to the toric-code
distance.

### Noise and success definitions

Each physical qubit independently receives an `X`, `Y`, or `Z` error with
probability `p/3` for each Pauli, for total physical error probability `p`.
The two decoders must use the same code convention and logical-failure
definition. A trial is a logical failure when the residual of the sampled
error and proposed correction has nontrivial logical action, not merely when
the correction has a syndrome inconsistent with the input.

### Logical-error sampling

The archived unitary-decouple curves used:

- `max_sim = 1_000_000_000` per physical error rate;
- `max_error = 2_000`;
- 120 worker processes;
- stopping when either limit is reached.

The historical distributed implementation assigns each worker a chunk error
limit of `ceil(remaining_error / worker_count)`. On the first round this is
17, so converged points finish with 2,040 failures. The release must preserve
this behavior for historical interpretation and document it in the result
metadata. A new deterministic runner may use a stricter global stop rule, but
its output must be labeled as a new run rather than an exact recreation of the
archived run.

The archived BP-OSD curve is assembled from single-point result files. For
each `(distance, p)` pair, the exporter selects the valid result with the
largest `nsim`, then sorts the selected points by `p`. The manifest records
every selected source filename and checksum so selection cannot silently
change.

### Confidence intervals

The figure uses the profile-likelihood binomial interval with likelihood ratio
parameter `h = 1000`, including the corresponding nonzero upper limit for
zero-failure observations. This method is implemented once in
DecoderBenchmarks and its normalized `low`, `estimate`, and `high` values are
exported to the companion dataset.

Wilson intervals implemented by other code in the companion are not used for
this historical figure. Documentation must name the method rather than call
the bars generic binomial error bars.

### Decoder configuration

The comparison intentionally preserves the historical asymmetry:

- BP-OSD is reconstructed for each point with a prior matched to that point's
  physical error rate;
- each toric-code sector of the Julia unitary-decouple decoder uses unweighted
  MWPM, with no edge
  probabilities passed to the matching graph.

This configuration is a historical fact, not a recommended fairness standard.
The figure caption and dataset metadata must state it. A future weighted-MWPM
comparison is a new experiment and is outside this release.

### Timing

Timing is reported as decoder-kernel time, not end-to-end benchmark latency.
For each point the historical runner:

1. compiled the decoder before entering the point loop;
2. executed 100 warm-up decodes;
3. executed 10,000 timed decodes;
4. timed only the call to `decode`;
5. excluded error sampling, syndrome generation, result validation, process
   startup, compilation, and file I/O.

The unitary-decouple timing was generated by the packed-bit Julia
implementation. It must not be attributed to the Sage/Python decoder exposed
elsewhere in the companion repository. The different physical-error grids
used for unitary-decouple and BP-OSD timing are retained and recorded
explicitly.

## Dependency and licensing boundary

TensorQEC and DecoderBenchmarks are MIT-licensed. The historical TensorQEC
branch currently refers to QECCore and SparseBlossom through local absolute
paths; those references are release blockers.

QECCore must be resolved through a public package version or a public URL and
immutable commit. SparseBlossom must become a separately installable public
dependency with:

- an explicit license for the Julia wrapper;
- an Apache-2.0 notice for the wrapped PyMatching code;
- a build process that does not rely on an untracked local shared library;
- a stable package UUID and version;
- clean-machine installation and test instructions;
- a pinned PyMatching source revision.

If the historical wrapper cannot be released with clear ownership and a
portable build, the logical-error portion may still be released with archived
data and a replacement matching backend, but the historical timing numbers
must then remain archived-only. Replacement-backend timing must be presented
as a new result and must not overwrite the article dataset.

## Public artifacts

### TensorQEC branch

The branch publishes:

- the minimal decoder implementation and tests;
- a small checked-in code fixture;
- a benchmark-independent example that compiles the decoder and corrects a
  deterministic low-weight error;
- a project environment with no path dependencies;
- machine-readable version output.

### DecoderBenchmarks branch

The branch publishes:

- `run_logical_error.jl` with smoke and full presets;
- `run_timing.jl` with the historical timing scope;
- the BP-OSD runner with explicit decoder parameters;
- deterministic seed derivation for process and point streams;
- the profile-likelihood implementation and tests;
- compact BB-code inputs for distances 4, 6, 8, and 10;
- the selected archived raw JSON files;
- `provenance.json` mapping raw files to normalized output;
- an exporter that creates the companion dataset without editing tracked
  source data in place.

Full generated output goes under ignored build directories. Only the curated
archived inputs, selected raw results, normalized export, and provenance
manifest are tracked.

### Companion repository

The companion publishes:

- an updated `decoding_benchmark.json` schema;
- a checksum manifest for the normalized dataset and external source refs;
- `make reproduce-decoder-smoke`;
- `make replot-decoder-benchmark`;
- release tests for metadata completeness, source hashes, confidence interval
  values, curve dimensions, and figure data selection;
- documentation separating archived results, deterministic smoke runs, and
  statistically reproducible full runs.

The normalized metadata includes:

- source repository, commit, and tag;
- data-generation and export dates as separate fields;
- Julia, Python, BP-OSD, TensorQEC, and matching-backend versions;
- OS, CPU, worker count, warm-up count, and timed sample count;
- code parameters and physical qubit counts;
- noise model, stop rule, seed status, interval method, and matching weights;
- timing inclusion and exclusion rules;
- selected raw source file hashes.

For the archived article run, `seed_status` is `not_recorded`; no replacement
seed is invented.

## Verification

Verification is layered so routine CI does not attempt the full Monte Carlo
run.

### Dependency gate

A clean temporary Julia depot must instantiate TensorQEC and run its focused
decoder tests without any local development packages or absolute paths.

### Deterministic smoke gate

The benchmark branch runs a fixed small matrix of distances, physical error
rates, and seeds. It verifies:

- syndrome consistency;
- logical success for selected correctable errors;
- stable result schema;
- deterministic counts for the fixed smoke configuration;
- identical profile-likelihood output for pinned test cases.

### Archived-data gate

The exporter verifies all source checksums, reproduces the normalized JSON,
and compares it byte-for-byte with the companion copy after canonical JSON
formatting. Tests pin representative logical-error and timing points,
including the 2,040-failure stop-rule behavior and the 5.317783355712891e-6
second distance-4 timing point.

### Companion release gate

The companion continues to run:

```bash
python3 -m unittest discover -s tests/release -v
python3 tools/verify_release.py
```

The smoke reproduction may download only public, immutable source refs. It
must write generated output under `build/reproduction/`.

## Rollout order

1. Make SparseBlossom independently licensed, buildable, and pinned.
2. Create and validate the TensorQEC paper branch and tag.
3. Create the DecoderBenchmarks paper branch from the historical commit.
4. Add deterministic benchmark entrypoints, curated inputs, archived raw
   results, statistics, and provenance export.
5. Tag the benchmark branch.
6. Update the companion dataset, documentation, Make targets, and tests to
   consume the two immutable tags.
7. Update the article caption and methods text to match the published
   provenance and timing contract.
8. Reply to issue #1 with exact tags, commands, limitations, and validation
   results.

## Deferred work

The following work is intentionally excluded from this release:

- merging the historical unitary-decouple decoder into the current TensorQEC
  main branch;
- weighted or correlated MWPM experiments;
- regenerating the article figure with a replacement decoder backend;
- performance comparisons across new hardware;
- opening the later bundle-BP matching experiments;
- claiming exact Monte Carlo replay for the unseeded historical run.

These items may be pursued after the historical release is stable, but they
must produce new datasets and must not change the archived article record.
