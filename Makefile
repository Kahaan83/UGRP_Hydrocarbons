# Minimal entry points (run with the venv python: make PYTHON=.venv/Scripts/python.exe test)
PYTHON ?= python

.PHONY: reproduce test
reproduce:   ## rebuild the frozen v5 scoreboard and check it against SHA256SUMS_v5_base.txt
	$(PYTHON) scripts/v5/base/reproduce_base.py

test:
	$(PYTHON) -m pytest -q
