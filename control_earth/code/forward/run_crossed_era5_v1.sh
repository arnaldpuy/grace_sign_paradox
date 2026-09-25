#!/bin/zsh
# Crossed-truth planet (D039): ERA5-Land background, GLDAS-Noah auxiliary.
# Builds the truth, forward-models it into synthetic Level-2 and runs the
# unchanged recovery pipeline on all three centres. ~10 h end to end.
set -e

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
PY="${PYTHON:-python3}"
SPEC=$ROOT/datasets/truth/truth_spec_crossed_era5_v1.json
PARAMS=$ROOT/datasets/truth/params_iter3.json

cd $ROOT
export CE_PLANET=crossed_era5_v1

echo "=== [1/3] truth builder (ERA5-Land background) $(date) ==="
$PY code/forward/01_truth_builder.py --spec $SPEC

echo "=== [2/3] forward operator $(date) ==="
$PY code/forward/02_forward_operator.py --spec $SPEC --params $PARAMS

echo "=== [3/3] recovery, three centres $(date) ==="
$PY code/forward/04_recover.py --centres CSR JPL GFZ

echo "=== CROSSED PLANET COMPLETE $(date) ==="
