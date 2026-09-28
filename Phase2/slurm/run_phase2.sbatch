#!/usr/bin/env bash
#SBATCH --mem=8G
#SBATCH --time=24:00:00

set -euo pipefail

stage=${1:?stage required}
list=${2:?input list required}
output=${3:?output directory required}
max_events=${4:--1}
index=$(( ${5:-0} + ${SLURM_ARRAY_TASK_ID:?submit as an array} ))

case "$stage" in
  rereco) cfg=run_rereco_cfg.py ;;
  ntuple) cfg=run_ntuple_cfg.py ;;
  *) echo "Unknown Phase2 stage: $stage" >&2; exit 2 ;;
esac

mapfile -t inputs < "$list"
input=${inputs[index]:-}
[[ -n $input ]] || { echo "Missing input for index $index" >&2; exit 1; }
case "$input" in
  root://*|file:/*) ;;
  /store/*) input="root://cms-xrd-global.cern.ch/$input" ;;
  /*) input="file:$input" ;;
  *) echo "Invalid input path: $input" >&2; exit 2 ;;
esac

printf -v number '%05d' "$index"
final="$output/${stage}_${number}.root"
partial="$output/.partial/${stage}_${number}_${SLURM_JOB_ID}.root"
if [[ -s $final ]]; then
  echo "Already exists: $final"
  exit 0
fi
mkdir -p "$output/.partial"
trap 'rm -f "$partial"' EXIT

proxy_args=()
if [[ $input == root://* ]]; then
  proxy="$HOME/.globus/cms-proxy"
  [[ -s $proxy ]] || { echo "Missing VOMS proxy: $proxy" >&2; exit 1; }
  proxy_args=(--bind "$proxy:$proxy:ro"
    --env "X509_USER_PROXY=$proxy"
    --env X509_CERT_DIR=/cvmfs/grid.cern.ch/etc/grid-security/certificates)
fi

base="$HOME/workspace/deepmuonreco/CMSSW_14_0_9"
[[ -f $base/.deepmuonreco-build-status && $(< "$base/.deepmuonreco-build-status") == ready ]] || {
  echo "CMSSW build is not ready: $base" >&2
  exit 1
}

echo "[$(date -Is)] START stage=$stage index=$index input=$input node=$(hostname)"
apptainer exec --cleanenv \
  --bind /cvmfs:/cvmfs \
  --bind "$HOME:$HOME" \
  "${proxy_args[@]}" \
  /cvmfs/unpacked.cern.ch/registry.hub.docker.com/cmssw/el8:x86_64 \
  bash -s -- "$base" "$cfg" "$input" "$partial" "$max_events" <<'IN_CONTAINER'
set -euo pipefail
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el8_amd64_gcc12
cd "$1/src"
eval "$(scram runtime -sh)"
cd "$(dirname "$4")"
cmsRun "$1/src/DeepMuonRecoSample/Phase2/test/$2" \
  "inputFiles=$3" "outputFile=$(basename "$4")" "maxEvents=$5"
IN_CONTAINER

test -s "$partial"
mv "$partial" "$final"
echo "[$(date -Is)] DONE $final"
