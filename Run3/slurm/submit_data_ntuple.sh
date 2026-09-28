#!/usr/bin/env bash
# Usage: bash Run3/slurm/submit_data_ntuple.sh [--test] [production-id-or-list.txt]
set -euo pipefail

here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd -- "$here/../.." && pwd)
base=/users/hep/joshin/workspace/deepmuonreco/CMSSW_14_0_21_patch1
image=/cvmfs/unpacked.cern.ch/registry.hub.docker.com/cmssw/el8:x86_64
output_root="$HOME/workspace/.store/deepmuonreco/chunk/Run3/Data"

mode=full
if [[ ${1:-} == --test ]]; then mode=test; shift; fi
if (( $# > 1 )); then
  echo "Usage: $0 [--test] [production-id-or-list.txt]" >&2
  exit 2
fi
selector=${1:-data-run3-muon0-2024cde-v001}
if [[ -f "$selector" ]]; then
  list=$(readlink -f "$selector")
elif [[ $selector == *.txt ]]; then
  list="$repo/Run3/inputs/$selector"
else
  list="$repo/Run3/inputs/$selector.txt"
fi
[[ -f $list ]] || { echo "Missing input list: $list" >&2; exit 2; }
[[ $list == *.txt ]] || { echo "Input list must be a .txt file" >&2; exit 2; }
tag=$(basename "$list" .txt)
[[ $tag =~ ^[a-z0-9][a-z0-9.-]*$ ]] || {
  echo "Input list name must use lowercase letters, numbers, dots or hyphens" >&2
  exit 2
}
[[ $mode == test ]] && tag=$tag-smoke
work="$output_root/$tag"
manifest="$work/jobs.tsv"
logs="$work/logs/ntuple"
output="$work/ntuple"
proxy=${X509_USER_PROXY:-$HOME/.globus/cms-proxy}
job_name="run3-data-ntuple-$tag"

command -v apptainer >/dev/null || { echo "Run this on a node with Apptainer" >&2; exit 1; }
[[ -f $base/.deepmuonreco-build-status && $(<"$base/.deepmuonreco-build-status") == ready ]] || {
  echo "Build CMSSW first: bash Run3/setup_cmssw.sh" >&2
  exit 1
}
[[ -r $proxy ]] || { echo "Missing VOMS proxy: $proxy" >&2; exit 1; }
voms-proxy-info -file "$proxy" -exists -valid 1:00 >/dev/null || {
  echo "Renew the CMS VOMS proxy before submitting" >&2
  exit 1
}

mkdir -p "$logs" "$output"
temporary=$(mktemp "$work/.jobs.XXXXXXXX")
trap 'rm -f "$temporary"' EXIT

apptainer exec --cleanenv \
  --bind /cvmfs:/cvmfs \
  --bind "$HOME:$HOME" \
  --env "X509_USER_PROXY=$proxy,X509_CERT_DIR=/cvmfs/grid.cern.ch/etc/grid-security/certificates,X509_VOMS_DIR=/cvmfs/grid.cern.ch/etc/grid-security/vomsdir" \
  "$image" bash -s -- "$repo" "$base" "$list" "$temporary" "$mode" <<'IN_CONTAINER'
set -euo pipefail
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el8_amd64_gcc12
cd "$2/src"
eval "$(scram runtime -sh)"
if [[ $5 == test ]]; then
  python3 "$1/Run3/slurm/plan_data_ntuple.py" --test "$3" > "$4"
else
  python3 "$1/Run3/slurm/plan_data_ntuple.py" "$3" > "$4"
fi
IN_CONTAINER

if [[ -e $manifest ]]; then
  cmp -s "$temporary" "$manifest" || {
    echo "Existing jobs.tsv differs; use a new input-list version" >&2
    exit 1
  }
else
  mv "$temporary" "$manifest"
fi

missing=()
index=0
selected=0
while IFS=$'\t' read -r events skip input filename; do
  [[ $events =~ ^[1-9][0-9]*$ && $skip =~ ^[0-9]+$ && $filename =~ ^ntuple-[0-9]{4,}\.root$ ]] || {
    echo "Invalid row $index in $manifest" >&2
    exit 1
  }
  if [[ -e $output/$filename ]]; then
    [[ -s $output/$filename ]] || { echo "Empty output: $output/$filename" >&2; exit 1; }
  else
    missing+=("$index")
  fi
  selected=$((selected + events))
  index=$((index + 1))
done < "$manifest"

echo "production: $tag"
echo "planned: $index jobs, $selected events"
echo "remaining: ${#missing[@]} jobs"
echo "output: $output"
if (( ${#missing[@]} == 0 )); then exit 0; fi
if [[ -n $(squeue -h -u "$USER" -n "$job_name" -o '%i') ]]; then
  echo "Jobs for this production are already queued" >&2
  exit 1
fi

for ((offset=0; offset<index; offset+=10000)); do
  batch=()
  for task in "${missing[@]}"; do
    (( task < offset )) && continue
    (( task >= offset + 10000 )) && break
    batch+=("$((task - offset))")
  done
  (( ${#batch[@]} > 0 )) || continue
  array=$(IFS=,; echo "${batch[*]}")
  sbatch --job-name="$job_name" --array="$array" \
    --output="$logs/%A_%a.out" \
    "$here/run_data_ntuple.sbatch" "$manifest" "$output" "$proxy" "$base" "$repo" "$offset"
done
