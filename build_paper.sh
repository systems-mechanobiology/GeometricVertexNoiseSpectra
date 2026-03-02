#!/usr/bin/env bash
set -euo pipefail

export TMPDIR="${TMPDIR:-$HOME/tmp}"
mkdir -p "$TMPDIR"

# Resolve repo root (directory containing this script)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

clean_aux() {
  rm -f paper.aux paper.bbl paper.blg paper.log paper.fls \
        paper.fdb_latexmk paper.out paper.toc paper.bcf paper.run.xml \
        paper.thm paper.nav paper.snm
}

MAX_ATTEMPTS=10

for attempt in $(seq 1 $MAX_ATTEMPTS); do
  clean_aux
  if quarto render paper.qmd --to pdf; then
    exit 0
  fi
  echo "Build attempt $attempt/$MAX_ATTEMPTS failed (TinyTeX may have auto-installed missing packages). Retrying..."
done

echo "Build failed after $MAX_ATTEMPTS attempts."
exit 1
