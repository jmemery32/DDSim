.PHONY: env install test example clean

# Create (or update) the conda environment and install ddsim into it, editable.
# Requires conda/miniconda/mamba on PATH.
env:
	conda env update -f environment.yml -n ddsim --prune
	@echo "Activate it with: conda activate ddsim"

# Install/reinstall ddsim (+ test dependencies) into the CURRENT Python environment,
# editable, without touching conda. Use this if you already have numpy/scipy/numba
# set up some other way.
install:
	pip install -e ".[dev]"

test:
	pytest

# Run the bundled example (a 20x20x20 cube, uniform stress) for a couple of nodes.
example:
	ddsim -base example1 -conpath examples/example1/ -parpath examples/example1/ \
		-v -doid_list 10,0 -scale 100

clean:
	find . -name '__pycache__' -not -path './reference/*' -exec rm -rf {} +
	rm -rf .pytest_cache build dist src/*.egg-info
