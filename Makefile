# Makefile for UGRP Hydrocarbons v2 Reproduction and Testing
PYTHON ?= python

.PHONY: reproduce test verify clean help v4

help:
	@echo "Available targets:"
	@echo "  make reproduce  - Run complete v2 reproduction pipeline (Phases A-E)"
	@echo "  make test       - Run all pytest test suites"
	@echo "  make verify     - Verify presence and integrity of all v2 deliverables"
	@echo "  make v4         - Run the v4 integrity reproduction (scripts/v4/reproduce_all.py, added in Task 19)"

reproduce:
	$(PYTHON) scripts/v2/reproduce_all.py

verify:
	$(PYTHON) scripts/v2/reproduce_all.py --verify-only

test:
	$(PYTHON) -m pytest tests/

v4:
	$(PYTHON) scripts/v4/reproduce_all.py
