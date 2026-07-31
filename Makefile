PYTHON ?= python3
SAGE ?= sage
JULIA ?= julia
DECODER_PYTHON ?= python3.13
DECODER_BENCHMARK_REPOSITORY ?= https://github.com/nzy1997/DecoderBenchmarks.jl.git
DECODER_BENCHMARK_REF ?= paper-decoder-benchmark-2026.2
DECODER_BENCHMARK_COMMIT ?= d3690ec0ab1018dcc4f49166715fecc5d3e6d5ab
override DECODER_BENCHMARK_DIR := build/reproduction/decoder-benchmark

.PHONY: check-julia-env test-python test-julia verify-manifest verify reproduce reproduce-full
.PHONY: prepare-decoder-benchmark verify-decoder-archive reproduce-decoder-smoke replot-decoder-benchmark

check-julia-env:
	$(JULIA) --startup-file=no --project=julia/ToricBuilder -e 'using Pkg; Pkg.instantiate(; allow_autoprecomp=false); using Oscar; println("Julia/Oscar environment: ok")'

test-python:
	cd python && $(SAGE) -pip install -e '.[test]'
	cd python && $(SAGE) -python -m compileall -q decoder_core decoupling unitary_decouple_based_decoder
	cd python && for test in tests/sage/decoupling/test_*.sage tests/sage/unitary/test_*.sage; do $(SAGE) "$$test" || exit 1; done
	cd python && $(SAGE) scripts/unitary_decouple_based_decoder/benchmark_paper_bb_family.sage --smoke
	cd python && $(SAGE) -c 'import decoder_core, decoupling, unitary_decouple_based_decoder'

test-julia:
	$(JULIA) --project=julia/ToricBuilder -e 'using Pkg; Pkg.instantiate(); Pkg.test()'

verify-manifest:
	$(PYTHON) tools/release_manifest.py --check

verify:
	$(PYTHON) -m unittest discover -s tests/release -v
	$(PYTHON) tools/verify_release.py
	$(MAKE) verify-manifest

reproduce:
	mkdir -p build/reproduction
	rm -rf build/reproduction/python
	mkdir -p build/reproduction/python
	cp -R python/decoder_core python/decoupling python/results python/scripts build/reproduction/python/
	cd build/reproduction/python && $(SAGE) scripts/decoupling/reproduce_666_color_code.sage
	cd build/reproduction/python && $(SAGE) scripts/decoupling/reproduce_488_color_code.sage
	cd build/reproduction/python && $(SAGE) scripts/decoupling/reproduce_bb_codes.sage --rows benchmark-01 --output-json ../python-bb-benchmark-01.json --output-csv ../python-bb-benchmark-01.csv --output-markdown ../python-bb-benchmark-01.md --checkpoint-json ../python-bb-benchmark-01.checkpoint.json
	$(JULIA) --project=julia/ToricBuilder julia/ToricBuilder/example/scripts/color_code.jl
	$(JULIA) --project=julia/ToricBuilder -e 'include("julia/ToricBuilder/example/scripts/decouple_bbcodes.jl"); main(list=ab_list1[1:1], results_path="build/reproduction/julia-bb-case-1.md", cache_dir="build/reproduction/julia-bb-cache", warmup=false, show_progress=false)'
	$(JULIA) --project=julia/ToricBuilder julia/ToricBuilder/example/scripts/plot_area_comparison.jl
	$(JULIA) --project=julia/ToricBuilder julia/ToricBuilder/example/scripts/plot_decoding_result_from_data.jl
	$(JULIA) --project=julia/ToricBuilder julia/ToricBuilder/example/scripts/transported_cnot_support.jl

reproduce-full:
	mkdir -p build/reproduction
	rm -rf build/reproduction/python
	mkdir -p build/reproduction/python
	cp -R python/decoder_core python/decoupling python/results python/scripts build/reproduction/python/
	cd build/reproduction/python && $(SAGE) scripts/decoupling/reproduce_bb_codes.sage --resume --output-json ../python-bb-full.json --output-csv ../python-bb-full.csv --output-markdown ../python-bb-full.md --checkpoint-json ../python-bb-full.checkpoint.json

prepare-decoder-benchmark:
	rm -rf $(DECODER_BENCHMARK_DIR)
	mkdir -p build/reproduction
	git clone --depth 1 --branch $(DECODER_BENCHMARK_REF) $(DECODER_BENCHMARK_REPOSITORY) $(DECODER_BENCHMARK_DIR)
	test "$$(git -C $(DECODER_BENCHMARK_DIR) rev-parse HEAD)" = "$(DECODER_BENCHMARK_COMMIT)"

verify-decoder-archive: prepare-decoder-benchmark
	cd $(DECODER_BENCHMARK_DIR) && $(JULIA) --project -e 'using Pkg; Pkg.instantiate(; allow_autoprecomp=false); Pkg.build("SparseBlossom")'
	cd $(DECODER_BENCHMARK_DIR) && $(JULIA) --project paper/export_normalized.jl $(abspath build/reproduction/decoder-export)
	cmp julia/ToricBuilder/results/decoding_benchmark.json build/reproduction/decoder-export/decoding_benchmark.json

reproduce-decoder-smoke: prepare-decoder-benchmark
	$(MAKE) -C $(DECODER_BENCHMARK_DIR) PYTHON=$(DECODER_PYTHON) paper-python-init
	$(MAKE) -C $(DECODER_BENCHMARK_DIR) paper-bposd-test
	$(MAKE) -C $(DECODER_BENCHMARK_DIR) JL='$(JULIA) --project' paper-smoke

replot-decoder-benchmark:
	$(JULIA) --project=julia/ToricBuilder julia/ToricBuilder/example/scripts/plot_decoding_result_from_data.jl
