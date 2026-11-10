#!/usr/bin/env bash
set -euo pipefail

# Run from the repository root. Use the same goose-doc checkout as docs.sh.
exec python3 "$(dirname "${GOOSE_DOC:-../goose-doc/main.goose}")/doctest.py" \
  --symbols api.txt --preamble tests/doc_preamble.goose --require-examples \
  "$@" radix_sort.goose iitree.goose lapper.goose fastx.goose
