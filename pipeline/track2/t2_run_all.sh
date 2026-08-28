#!/usr/bin/env bash
# ==============================================================================
# t2_run_all.sh - Track 2 drug repositioning pipeline, end to end
#
#   bash pipeline/track2/t2_run_all.sh
#
# Track 2 starts FROM the Track 1 result by design: the task is to go from
# variant and mechanism to candidate medicines. The causal gene is an input
# here, not a leaked answer.
# ==============================================================================
set -euo pipefail
D="$(cd "$(dirname "$0")" && pwd)"
source "$D/../00_config.sh"

log "===== T2-01 target network ====="
python3 "$D/t2_01_target_network.py" > /dev/null

log "===== T2-02 drug-gene evidence ====="
python3 "$D/t2_02_drug_evidence.py" > /dev/null

log "===== T2-03 mechanistic filter ====="
python3 "$D/t2_03_mechanism_filter.py" > /dev/null

log "===== T2-04 readthrough branch ====="
python3 "$D/t2_04_readthrough_branch.py" > /dev/null

log "===== mirroring results ====="
bash "$D/../sync_results.sh" > /dev/null

log "TRACK 2 PIPELINE COMPLETE"
echo
echo "  results in \$RESULTS/t2_0*.txt"
