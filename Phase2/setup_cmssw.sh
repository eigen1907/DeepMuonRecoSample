#!/usr/bin/env bash
# Run once in a shell that has Apptainer (an interactive compute node here).
set -euo pipefail

if ! command -v apptainer >/dev/null 2>&1; then
  echo "Apptainer is unavailable on $(hostname). Use an interactive compute shell." >&2
  exit 1
fi

repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
base="$HOME/workspace/deepmuonreco/CMSSW_14_0_9"

apptainer exec --cleanenv \
  --bind /cvmfs:/cvmfs \
  --bind "$HOME:$HOME" \
  /cvmfs/unpacked.cern.ch/registry.hub.docker.com/cmssw/el8:x86_64 \
  bash -s -- "$repo" "$base" <<'IN_CONTAINER'
set -euo pipefail
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el8_amd64_gcc12

repo=$1
base=$2
cd "$(dirname "$base")"
if [[ ! -d "$base/src" ]]; then
  scram project CMSSW CMSSW_14_0_9
fi
cd "$base/src"
if [[ ! -e DeepMuonRecoSample ]]; then
  ln -s "$repo" DeepMuonRecoSample
elif [[ $(realpath DeepMuonRecoSample) != "$repo" ]]; then
  echo "DeepMuonRecoSample points to a different repository" >&2
  exit 1
fi
status="$base/.deepmuonreco-build-status"
printf 'building\n' > "$status"
eval "$(scram runtime -sh)"
scram b -j 4
printf 'ready\n' > "$status"
IN_CONTAINER

echo "DONE $base"
