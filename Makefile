# Makefile for UGRP Hydrocarbons v2 Reproduction and Testing
PYTHON ?= python

.PHONY: reproduce test verify clean help

help:
	@echo "Available targets:"
	@echo "  make reproduce  - Run complete v2 reproduction pipeline (Phases A-E)"
	@echo "  make test       - Run all pytest test suites"
	@echo "  make verify     - Verify presence and integrity of all v2 deliverables"

reproduce:
	$(PYTHON) scripts/v2/reproduce_all.py

verify:
	$(PYTHON) scripts/v2/reproduce_all.py --verify-only

test:
	$(PYTHON) -m pytest tests/
