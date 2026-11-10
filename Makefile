.PHONY: help test test-full lint format check smoke

# pytest-xdist worker count; set PYTEST_WORKERS=1 to run serially when reading a traceback.
PYTEST_WORKERS ?= 3
# Extra pytest arguments, e.g. make test PYTEST_ARGS="-k ledger".
PYTEST_ARGS ?=
# One BLAS thread per worker; the workers already use the cores.
export OPENBLAS_NUM_THREADS := 1
export OMP_NUM_THREADS := 1
export MKL_NUM_THREADS := 1
export MPLBACKEND := Agg

help:
	@echo "test (fast tier, <120 s)  test-full (<300 s)  lint  format  check (lint + full)  smoke"

# Fast tier: under 120 s of test-process wall time, no single test over about 20 s.
test:
	python scripts/run_tests.py fast -n $(PYTEST_WORKERS) $(PYTEST_ARGS)

# Full tier: everything, under 300 s, no single slow test over about 90 s. The gate.
test-full:
	python scripts/run_tests.py full -n $(PYTEST_WORKERS) $(PYTEST_ARGS)

lint:
	python -m black --check --diff src tests scripts

format:
	python -m black src tests scripts

check: lint test-full

# Clean-environment install and import, outside the source tree.
smoke:
	python scripts/smoke.py
