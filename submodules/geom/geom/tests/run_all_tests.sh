#!/bin/bash

# Get the absolute path of the script directory (tests/)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Navigate to the tests directory
cd "$SCRIPT_DIR" || exit 1

# Run pytest
pytest -v --tb=short --durations=0

