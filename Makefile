PY := python
PIP := pip

help:
	@echo "Commands:"
	@echo ""
	@echo "    pep8          check for PEP8 style compliance"
	@echo ""

pep8:
	pep8 --show-source --ignore=W503,E226,E241 Conversion.py
