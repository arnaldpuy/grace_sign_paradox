#!/bin/zsh
# Recovery only, for the crossed-truth planet (D045).
#
# The truth build and the forward operator are both idempotent and already
# complete (datasets/truth/truth_tws_crossed_era5_v1.nc,
# datasets/synthetic_l2/crossed_era5_v1/{csr,gfz,jpl,auxiliary}), so a
# restart after an interrupted run needs this step alone. 04_recover.py takes
# --centres, so a centre that already produced
# datasets/output/tier_a_synth_crossed_era5_v1/<centre>/ can be dropped from
# the list instead of being recomputed.
#
# Launch detached so that it survives the parent shell:
#   nohup setsid zsh code/forward/run_crossed_recovery.sh > LOG 2>&1 &
set -e

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
PY="${PYTHON:-python3}"

cd $ROOT
export CE_PLANET=crossed_era5_v1

if [ $# -eq 0 ]; then
  set -- CSR JPL GFZ
fi

echo "=== recovery, centres: $@  $(date) ==="
caffeinate -i $PY code/forward/04_recover.py --centres "$@"
echo "=== CROSSED PLANET RECOVERY COMPLETE $(date) ==="
