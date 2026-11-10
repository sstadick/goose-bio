#!/usr/bin/env bash
set -euo pipefail

# Run from the repository root. GOOSE_DOC names goose-doc's main.goose.
exec "${GOOSE:-goose}" run "${GOOSE_DOC:-../goose-doc/main.goose}" -- \
  --title "goose-bio API" --symbols api.txt --require-docs --require-examples \
  --output "${DOC_OUTPUT:-API.md}" "$@" radix_sort.goose iitree.goose lapper.goose fastx.goose
