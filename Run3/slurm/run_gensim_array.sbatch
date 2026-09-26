#!/usr/bin/env bash
#SBATCH --mem=8G
#SBATCH --time=24:00:00
# Exclude nodes with broken CVMFS access as of 2026-09-26.
#SBATCH --exclude=entei,ho-oh

set -euo pipefail
stage=${1:?gensim or minbias required}
events=${2:?events per job required}
output=${3:?output directory required}
index=$(( ${4:-0} + ${SLURM_ARRAY_TASK_ID:?submit as an array} ))
case "$stage" in
  gensim) cfg=run_mc_signal_gensim_cfg.py ;;
  minbias) cfg=run_mc_minbias_gensim_cfg.py ;;
  *) echo "Unknown stage: $stage" >&2; exit 2 ;;
esac
printf -v number '%05d' "$index"
final="$output/${stage}_$number.root"
partial="$output/.partial/${stage}_${number}_${SLURM_JOB_ID}.root"

if [[ -e $final ]]; then
  echo "Already exists: $final"
  exit 0
fi
mkdir -p "$output/.partial"
echo "[$(date -Is)] START stage=$stage index=$index events=$events node=$(hostname)"

repo="$HOME/workspace/deepmuonreco/DeepMuonRecoSample"
apptainer exec --cleanenv \
  --bind /cvmfs:/cvmfs \
  --bind "$HOME/workspace:$HOME/workspace" \
  /cvmfs/unpacked.cern.ch/registry.hub.docker.com/cmssw/el8:x86_64 \
  bash -s -- "$repo" "$cfg" "$events" "$index" "$partial" <<'IN_CONTAINER'
set -euo pipefail
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el8_amd64_gcc12
cd /cvmfs/cms.cern.ch/el8_amd64_gcc12/cms/cmssw-patch/CMSSW_14_0_21_patch1/src
eval "$(scram runtime -sh)"
cd "$(dirname "$5")"
cmsRun "$1/Run3/test/$2" \
  maxEvents="$3" jobIndex="$4" outputFile="$(basename "$5")"
IN_CONTAINER

test -s "$partial"
if [[ -e $final ]]; then
  echo "Output appeared while running: $final" >&2
  exit 1
fi
mv "$partial" "$final"
echo "[$(date -Is)] DONE $final"
