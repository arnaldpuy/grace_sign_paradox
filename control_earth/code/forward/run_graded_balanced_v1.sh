#!/bin/zsh
# Sign-balanced graded planet (D048): GLDAS background, re-randomised signs.
# Truth -> forward -> recovery on all three centres. ~15 h end to end.
set -e

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
PY="${PYTHON:-python3}"
SPEC=$ROOT/datasets/truth/truth_spec_graded_balanced_v1.json
PARAMS=$ROOT/datasets/truth/params_iter3.json

cd $ROOT
export CE_PLANET=graded_balanced_v1

echo "=== [1/3] truth builder (GLDAS background) $(date) ==="
$PY code/forward/01_truth_builder.py --spec $SPEC

echo "=== [2/3] forward operator $(date) ==="
$PY code/forward/02_forward_operator.py --spec $SPEC --params $PARAMS

echo "=== [3/3] recovery, three centres $(date) ==="
$PY code/forward/04_recover.py --centres CSR JPL GFZ

echo "=== BALANCED PLANET COMPLETE $(date) ==="
