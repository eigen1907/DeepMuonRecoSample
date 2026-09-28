#!/usr/bin/env bash
# Run once on a node with Apptainer to build the local Run 3 Ntuplizer plugin.
set -euo pipefail

repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
base=/users/hep/joshin/workspace/deepmuonreco/CMSSW_14_0_21_patch1
image=/cvmfs/unpacked.cern.ch/registry.hub.docker.com/cmssw/el8:x86_64

command -v apptainer >/dev/null || { echo "Apptainer is required on this node" >&2; exit 1; }
mkdir -p "$(dirname "$base")"

apptainer exec --cleanenv \
  --bind /cvmfs:/cvmfs \
  --bind "$HOME:$HOME" \
  "$image" bash -s -- "$repo" "$base" <<'IN_CONTAINER'
set -euo pipefail
repo=$1
base=$2

source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el8_amd64_gcc12
if [[ ! -d "$base/src" ]]; then
  cd "$(dirname "$base")"
  scram project CMSSW CMSSW_14_0_21_patch1
fi

package="$base/src/DeepMuonRecoSample"
if [[ -e "$package" || -L "$package" ]]; then
  [[ $(readlink -f "$package") == "$repo" ]] || {
    echo "Different package already exists at $package" >&2
    exit 1
  }
else
  ln -s "$repo" "$package"
fi

cd "$base/src"
eval "$(scram runtime -sh)"
printf 'building\n' > "$base/.deepmuonreco-build-status"
scram b -j 4
printf 'ready\n' > "$base/.deepmuonreco-build-status"
echo "DONE: built $package"
IN_CONTAINER
