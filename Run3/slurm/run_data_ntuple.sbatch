#!/usr/bin/env bash
#SBATCH --mem=4G
#SBATCH --time=08:00:00
set -euo pipefail

manifest=${1:?jobs.tsv required}
output=${2:?output directory required}
proxy=${3:?VOMS proxy required}
base=${4:?CMSSW base required}
repo=${5:?repository path required}
task=$(( ${6:-0} + ${SLURM_ARRAY_TASK_ID:?submit as a Slurm array} ))
[[ $task =~ ^[0-9]+$ ]] || { echo "Invalid array index: $task" >&2; exit 2; }
mapfile -t rows < "$manifest"
(( task < ${#rows[@]} )) || { echo "Array index out of range: $task" >&2; exit 2; }
IFS=$'\t' read -r events skip input filename <<< "${rows[$task]}"
[[ $events =~ ^[1-9][0-9]*$ && $skip =~ ^[0-9]+$ && $filename =~ ^ntuple-[0-9]{4,}\.root$ ]] || {
  echo "Invalid jobs.tsv row for index $task" >&2
  exit 2
}

final="$output/$filename"
if [[ -e $final ]]; then
  [[ -s $final ]] || { echo "Empty output: $final" >&2; exit 1; }
  echo "Already complete: $final"
  exit 0
fi
[[ -r $proxy ]] || { echo "Missing VOMS proxy: $proxy" >&2; exit 1; }
[[ -f $base/.deepmuonreco-build-status && $(<"$base/.deepmuonreco-build-status") == ready ]] || {
  echo "CMSSW Ntuplizer is not built: $base" >&2
  exit 1
}

mkdir -p "$output/.partial"
partial="$output/.partial/${filename%.root}_${SLURM_JOB_ID:?}.root"
trap 'rm -f "$partial"' EXIT
echo "[$(date -Is)] START task=$task events=$events skip=$skip node=$(hostname)"

apptainer exec --cleanenv \
  --bind /cvmfs:/cvmfs \
  --bind "$HOME:$HOME" \
  --env "X509_USER_PROXY=$proxy,X509_CERT_DIR=/cvmfs/grid.cern.ch/etc/grid-security/certificates,X509_VOMS_DIR=/cvmfs/grid.cern.ch/etc/grid-security/vomsdir" \
  /cvmfs/unpacked.cern.ch/registry.hub.docker.com/cmssw/el8:x86_64 \
  bash -s -- "$repo" "$base" "$input" "$events" "$skip" "$partial" <<'IN_CONTAINER'
set -euo pipefail
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el8_amd64_gcc12
cd "$2/src"
eval "$(scram runtime -sh)"
cd "$(dirname "$6")"
cmsRun "$1/Run3/test/run_ntuple_cfg.py" \
  inputFiles="$3" maxEvents="$4" skipEvents="$5" \
  isMC=false outputFile="$6"
IN_CONTAINER

[[ -s $partial ]] || { echo "No ntuple produced: $partial" >&2; exit 1; }
[[ ! -e $final ]] || { echo "Output appeared while running: $final" >&2; exit 1; }
mv "$partial" "$final"
echo "[$(date -Is)] DONE $final"
