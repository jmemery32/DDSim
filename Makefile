.PHONY: venv install test example docs clean

# Create (or reuse) a plain venv at .venv/ and install ddsim into it, editable,
# with every optional extra (dev/exodus/plotting). No conda needed -- pure pip.
# Uses `python3` on PATH by default; override with e.g. `make venv PY=python3.11`
# if that's too old (>=3.10 required -- a plain `python3` is often an older
# system/distro Python; get a newer one via pyenv, Homebrew, or python.org,
# none of which is conda).
PY ?= python3
venv:
	@$(PY) -c 'import sys; sys.exit(0 if sys.version_info >= (3, 10) else 1)' || \
		{ echo "$(PY) is $$($(PY) --version 2>&1), but ddsim needs >=3.10."; \
		  echo "Point PY at a newer interpreter, e.g.: make venv PY=python3.11"; \
		  exit 1; }
	$(PY) -m venv .venv
	./.venv/bin/pip install --upgrade pip
	./.venv/bin/pip install -e ".[dev,exodus,plotting]"
	@echo "Activate it with: source .venv/bin/activate"

# Install/reinstall ddsim (+ test dependencies) into the CURRENT Python
# environment, editable. Use this instead of `make venv` if you already have
# some other Python environment (conda or otherwise) you'd rather use.
install:
	pip install -e ".[dev]"

test:
	pytest

# Run the bundled example (a 20x20x20 cube, uniform stress) for a couple of nodes.
example:
	ddsim -base example1 -conpath examples/example1/ -parpath examples/example1/ \
		-v -doid_list 10,0 -scale 100

docs:
	sphinx-build -b html docs/source docs/build/html
	@echo "Open docs/build/html/index.html"

clean:
	find . -name '__pycache__' -not -path './reference/*' -exec rm -rf {} +
	rm -rf .pytest_cache build dist src/*.egg-info docs/build docs/source/api/generated
