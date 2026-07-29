PYTHON ?= python3
SAGE ?= sage
JULIA ?= julia

.PHONY: test-python test-julia verify-manifest verify reproduce reproduce-full

test-python:
	cd python && $(SAGE) -pip install -e '.[test]'
	cd python && for test_file in tests/python/test_*.py; do $(SAGE) -python "$$test_file" || exit 1; done
	cd python/tests/sage/core && for test_file in test_*.sage; do $(SAGE) "$$test_file" || exit 1; done
	cd python/tests/sage/iso && for test_file in test_*.sage; do $(SAGE) "$$test_file" || exit 1; done
	cd python/tests/sage/unitary && for test_file in test_*.sage; do $(SAGE) "$$test_file" || exit 1; done

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
	cd python && $(SAGE) -python -m jupyter nbconvert --to notebook --execute results/isomorphism/666_color_code.ipynb --output-dir ../build/reproduction --output 666_color_code.executed.ipynb
	cd python && $(SAGE) -python -m jupyter nbconvert --to notebook --execute results/isomorphism/488_color_code.ipynb --output-dir ../build/reproduction --output 488_color_code.executed.ipynb
	cd python && $(SAGE) scripts/iso/validate_bb_instances.sage --rows 1 --qca-check skip --row-timeout-seconds 300 --stop-on-failure --output-json ../build/reproduction/python-bb-row-1.json --output-md ../build/reproduction/python-bb-row-1.md
	$(JULIA) --project=julia/ToricBuilder julia/ToricBuilder/example/scripts/color_code.jl
	$(JULIA) --project=julia/ToricBuilder -e 'include("julia/ToricBuilder/example/scripts/decouple_bbcodes.jl"); main(list=ab_list1[1:1], results_path="build/reproduction/julia-bb-case-1.md", cache_dir="build/reproduction/julia-bb-cache", warmup=false, show_progress=false)'
	$(JULIA) --project=julia/ToricBuilder julia/ToricBuilder/example/scripts/plot_area_comparison.jl
	$(JULIA) --project=julia/ToricBuilder julia/ToricBuilder/example/scripts/plot_decoding_result_from_data.jl
	$(JULIA) --project=julia/ToricBuilder julia/ToricBuilder/example/scripts/transported_cnot_support.jl

reproduce-full:
	mkdir -p build/reproduction
	cd python && $(SAGE) scripts/iso/validate_bb_instances.sage --qca-check product --row-timeout-seconds 1800 --stop-on-failure --output-json ../build/reproduction/python-bb-full.json --output-md ../build/reproduction/python-bb-full.md
