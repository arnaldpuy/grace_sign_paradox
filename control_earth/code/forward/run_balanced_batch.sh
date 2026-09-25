#!/bin/zsh
# Balanced layouts 2, 3 and the sign-reversed counterpart of layout 1
# (AMENDMENT_01_sign_balance.md section 3; D051). Sequential: truth ->
# forward -> recovery on all three centres per planet, ~14.5 h each.
set -e

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
PY="${PYTHON:-python3}"
PARAMS=$ROOT/datasets/truth/params_iter3.json

cd $ROOT

for PLANET in graded_balanced_v2 graded_balanced_v3 graded_balanced_v1_reversed; do
  SPEC=$ROOT/datasets/truth/truth_spec_$PLANET.json
  export CE_PLANET=$PLANET

  echo "=== [$PLANET 1/3] truth builder $(date) ==="
  $PY code/forward/01_truth_builder.py --spec $SPEC

  echo "=== [$PLANET 2/3] forward operator $(date) ==="
  $PY code/forward/02_forward_operator.py --spec $SPEC --params $PARAMS

  echo "=== [$PLANET 3/3] recovery, three centres $(date) ==="
  $PY code/forward/04_recover.py --centres CSR JPL GFZ

  echo "=== $PLANET COMPLETE $(date) ==="
done

echo "=== BALANCED BATCH COMPLETE $(date) ==="
