SHELL := /bin/bash

PY := python
PIP := pip

ENV_NAME := h2-kmc

.PHONY: help setup install test run plot lint docker-build docker-run clean

help:
	@echo "Targets:"
	@echo "  setup         Create conda env from environment.yml"
	@echo "  install       Ensure pip tools installed in current env"
	@echo "  test          Run unit tests"
	@echo "  run           Run parameter sweep (uses config.yaml)"
	@echo "  plot          Generate plots from results CSV"
	@echo "  lint          Basic lint: check syntax"
	@echo "  docker-build  Build Docker image"
	@echo "  docker-run    Run in Docker container"
	@echo "  clean         Remove results and plots"

setup:
	@echo "Create env: $(ENV_NAME)"
	@conda env create -f environment.yml || conda env update -f environment.yml
	@echo "Activate with: conda activate $(ENV_NAME)"

install:
	$(PIP) --version >/dev/null

test:
	$(PY) -m pytest -q

run:
	$(PY) run_sweep.py

plot:
	$(PY) analysis_and_plotting.py

lint:
	$(PY) -m py_compile kmc_simulation.py run_sweep.py analysis_and_plotting.py scientific_data.py || true

docker-build:
	docker build -t h2-kmc:latest .

docker-run:
	docker run --rm -v $(PWD):/app h2-kmc:latest

clean:
	rm -f results/*.csv || true
	rm -rf results/plots || true