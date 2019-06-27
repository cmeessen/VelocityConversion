PY := python
PIP := pip

help:
	@echo "Commands:"
	@echo ""
	@echo "    pycodestyle          check for code style conventions"
	@echo ""

pycodestyle:
	pycodestyle --show-source --show-pep8 --ignore=W503,E226,E241 VelocityConversion/__init__.py

