#!/usr/bin/env bash
#SBATCH --mem=8G
#SBATCH --time=24:00:00

set -euo pipefail
stage=${1:?processing stage required}
events=${2:?event count required}
output=${3:?output directory required}
index=$(( ${4:?array offset required} + ${SLURM_ARRAY_TASK_ID:?submit as an array} ))
input_dir=${5:?input directory required}
secondary_dir=${6:-}
ratio=${7:-0}
repo=${8:?repository path required}

case $stage in
  premix) source_stage=gensim; cfg=run_mc_premix_cfg.py ;;
  digiraw) source_stage=gensim; cfg=run_mc_digiraw_premix_cfg.py ;;
  reco) source_stage=digiraw; cfg=run_mc_reco_cfg.py ;;
  ntuple) source_stage=reco; cfg=run_ntuple_cfg.py ;;
  *) echo "Unknown processing stage: $stage" >&2; exit 2 ;;
esac
printf -v number '%05d' "$index"
input="$input_dir/${source_stage}_$number.root"
final="$output/${stage}_$number.root"
partial="$output/.partial/${stage}_${number}_${SLURM_JOB_ID}.root"
[[ -s $input ]] || { echo "Missing or empty input: $input" >&2; exit 1; }
if [[ -s $final ]]; then
  echo "Already exists: $final"
  exit 0
fi
[[ ! -e $final ]] || { echo "Empty output exists: $final" >&2; exit 1; }
mkdir -p "$output/.partial"
secondary=
cleanup() {
  if [[ -e $partial ]]; then rm -f -- "$partial"; fi
  if [[ $stage == premix && -n $secondary && -e $secondary ]]; then
    rm -f -- "$secondary"
  fi
  return 0
}
trap cleanup EXIT

if [[ $stage == premix ]]; then
  [[ $ratio =~ ^[1-9][0-9]*$ ]] || { echo "Invalid MinBias ratio: $ratio" >&2; exit 2; }
  secondary="$output/.partial/minbias_${number}_${SLURM_JOB_ID}.txt"
  for ((j=0; j<ratio; j++)); do
    mb_index=$(( index * ratio + j ))
    printf -v mb_number '%05d' "$mb_index"
    mb_file="$secondary_dir/minbias_$mb_number.root"
    [[ -s $mb_file ]] || { echo "Missing or empty MinBias input: $mb_file" >&2; exit 1; }
    printf 'file:%s\n' "$mb_file"
  done > "$secondary"
elif [[ $stage == digiraw ]]; then
  secondary="$secondary_dir/premix_$number.root"
  [[ -s $secondary ]] || { echo "Missing or empty premix input: $secondary" >&2; exit 1; }
fi

if [[ $stage == ntuple ]]; then
  cmssw_base="$HOME/workspace/deepmuonreco/CMSSW_14_0_21_patch1"
  status="$cmssw_base/.deepmuonreco-build-status"
  build_status=
  if [[ -f $status ]]; then read -r build_status < "$status"; fi
  [[ $build_status == ready ]] || {
    echo "Run 3 CMSSW plugin is not built. Run bash Run3/setup_cmssw.sh once." >&2
    exit 1
  }
else
  cmssw_base=/cvmfs/cms.cern.ch/el8_amd64_gcc12/cms/cmssw-patch/CMSSW_14_0_21_patch1
fi

echo "[$(date -Is)] START stage=$stage index=$index events=$events input=$input node=$(hostname)"
apptainer exec --cleanenv \
  --bind /cvmfs:/cvmfs \
  --bind "$HOME:$HOME" \
  /cvmfs/unpacked.cern.ch/registry.hub.docker.com/cmssw/el8:x86_64 \
  bash -s -- "$repo" "$cfg" "$stage" "$events" "$index" "$input" "$partial" "$secondary" "$cmssw_base" <<'IN_CONTAINER'
set -euo pipefail
repo=$1
cfg=$2
stage=$3
events=$4
index=$5
input=$6
partial=$7
secondary=$8
cmssw_base=$9
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el8_amd64_gcc12
cd "$cmssw_base/src"
eval "$(scram runtime -sh)"
cd "$(dirname "$partial")"
case $stage in
  premix)
    cmsRun "$repo/Run3/test/$cfg" minbiasList="$secondary" \
      maxEvents="$events" jobIndex="$index" outputFile="$partial"
    ;;
  digiraw)
    cmsRun "$repo/Run3/test/$cfg" inputFiles="file:$input" \
      premixFile="file:$secondary" maxEvents="$events" \
      jobIndex="$index" outputFile="$partial"
    ;;
  reco)
    cmsRun "$repo/Run3/test/$cfg" inputFiles="file:$input" \
      maxEvents="$events" outputFile="$partial"
    ;;
  ntuple)
    cmsRun "$repo/Run3/test/$cfg" inputFiles="file:$input" \
      maxEvents="$events" isMC=true isPremix=true outputFile="$partial"
    ;;
esac
IN_CONTAINER

test -s "$partial"
if [[ -e $final ]]; then
  echo "Output appeared while running: $final" >&2
  exit 1
fi
mv "$partial" "$final"
echo "[$(date -Is)] DONE $final"
