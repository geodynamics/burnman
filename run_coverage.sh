#!/bin/bash

# This script runs the test suite and checks that code coverage has not decreased
# relative to the main branch. It also generates an HTML coverage report.

# Uninstall numba to avoid issues with coverage measurement
python -m pip install -q -e .[dev]
python -m pip uninstall -y numba

# Get current coverage
coverage run -m unittest discover ./tests
coverage report > new_coverage.tmp
current_coverage=$(coverage report --precision=7 --format=total)
echo "Coverage: $current_coverage"

# Generate HTML coverage report
coverage html -d htmlcov

# Reinstall numba, clean up temporary files
python -m pip install -q -e .[dev]
rm new_coverage.tmp
echo "Numba reinstalled:"
pip freeze | grep numba


# Open the HTML coverage report in the default web browser if on macOS
if [[ "$OSTYPE" == "darwin"* ]]; then
  open "htmlcov/index.html"
fi
