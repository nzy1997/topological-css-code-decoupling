PYTHON ?= python3
SAGE ?= sage
JULIA ?= julia

.PHONY: test-python test-julia verify-manifest verify reproduce reproduce-full

test-python:
	cd python && $(SAGE) -pip install -e '.[test]'
	cd python && $(SAGE) -python -m compileall -q decoder_core decoupling unitary_decouple_based_decoder
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
