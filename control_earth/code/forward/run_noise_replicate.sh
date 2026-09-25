#!/bin/zsh
# Noise replicate of graded_balanced_v3 (amendment 01 section 5; D055/D056).
# Waits for the running balanced batch to finish, then truth -> forward ->
# three-centre recovery for graded_balanced_v3_noise2 (~14.5 h).
set -e

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
PY="${PYTHON:-python3}"
PLANET=graded_balanced_v3_noise2
SPEC=$ROOT/datasets/truth/truth_spec_$PLANET.json
PARAMS=$ROOT/datasets/truth/params_iter3.json

cd $ROOT
until grep -q "BALANCED BATCH COMPLETE" logs/recover_balanced_batch.log; do sleep 300; done
export CE_PLANET=$PLANET

echo "=== [$PLANET 1/3] truth builder $(date) ==="
$PY code/forward/01_truth_builder.py --spec $SPEC

echo "=== [$PLANET 2/3] forward operator $(date) ==="
$PY code/forward/02_forward_operator.py --spec $SPEC --params $PARAMS

echo "=== [$PLANET 3/3] recovery, three centres $(date) ==="
$PY code/forward/04_recover.py --centres CSR JPL GFZ

echo "=== NOISE REPLICATE COMPLETE $(date) ==="
