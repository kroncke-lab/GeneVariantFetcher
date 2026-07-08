# GeneVariantFetcher — common entry points for collaborators.
# Run `make help` for the list. See docs/QUICKSTART.md for the full guide.

.DEFAULT_GOAL := help
PY := .venv/bin/python
EMAIL ?= brett.kroncke@gmail.com
OUTPUT ?= ./results

.PHONY: help install test run

help:  ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) \
		| awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-10s\033[0m %s\n", $$1, $$2}'

install:  ## Create .venv and install GVF with browser + dev extras
	python3.11 -m venv .venv
	$(PY) -m pip install -e ".[browser,dev]"
	$(PY) -m playwright install chromium

test:  ## Run the offline unit suite (the CI-gated set)
	$(PY) -m pytest tests/unit -q

run:  ## Run the pipeline on a gene: make run GENE=KCNH2 [EMAIL=brett.kroncke@gmail.com] [OUTPUT=./results]
	@test -n "$(GENE)" || { echo "Usage: make run GENE=<SYMBOL> [EMAIL=brett.kroncke@gmail.com] [OUTPUT=./results]"; exit 2; }
	$(PY) -m cli gvf-run $(GENE) --email $(EMAIL) --output $(OUTPUT)
