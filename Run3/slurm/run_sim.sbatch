#!/usr/bin/env bash
#SBATCH --job-name=run3-sim
#SBATCH --mem=8G
#SBATCH --time=08:00:00
# Exclude nodes with broken CVMFS access as of 2026-09-26.
#SBATCH --exclude=entei,ho-oh
#SBATCH --output=slurm-%x-%j.out

set -euo pipefail

mode=${1:-gensim}
case "$mode" in
  gensim|reco) ;;
  *) echo "Usage: sbatch $0 [gensim|reco]" >&2; exit 2 ;;
esac

repo="$HOME/workspace/deepmuonreco/DeepMuonRecoSample"
out="$HOME/workspace/deepmuonreco/run3-slurm-output/$SLURM_JOB_ID"
image="/cvmfs/unpacked.cern.ch/registry.hub.docker.com/cmssw/el8:x86_64"
mkdir -p "$out"

apptainer exec --cleanenv \
  --bind /cvmfs:/cvmfs \
  --bind "$HOME/workspace:$HOME/workspace" \
  "$image" bash -s -- "$out" "$repo" "$mode" <<'IN_CONTAINER'
set -euo pipefail
out=$1
repo=$2
mode=$3

source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el8_amd64_gcc12
cd /cvmfs/cms.cern.ch/el8_amd64_gcc12/cms/cmssw-patch/CMSSW_14_0_21_patch1/src
eval "$(scram runtime -sh)"
cd "$out"

cmsRun "$repo/Run3/test/run_mc_signal_gensim_cfg.py" \
  maxEvents=1 jobIndex=0 outputFile=gensim.root

if [[ "$mode" == reco ]]; then
  cmsRun "$repo/Run3/test/run_mc_minbias_gensim_cfg.py" \
    maxEvents=1000 jobIndex=0 outputFile=minbias.root
  cmsRun "$repo/Run3/test/run_mc_digiraw_cfg.py" \
    inputFiles=file:gensim.root secondaryInputFiles=file:minbias.root \
    maxEvents=1 jobIndex=0 outputFile=digiraw.root
  cmsRun "$repo/Run3/test/run_mc_reco_cfg.py" \
    inputFiles=file:digiraw.root maxEvents=1 outputFile=reco.root
fi
IN_CONTAINER

ls -lh "$out"
