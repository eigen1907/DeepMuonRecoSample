#!/usr/bin/env bash
# Submit one Run 3 MC stage at a time; GEN-SIM and MinBias use submit_gensim.sh.
# Usage: bash submit_processing.sh premix TAG GENSIM_DIR MINBIAS_DIR [--test]
#        bash submit_processing.sh digiraw TAG GENSIM_DIR [--test]
#        bash submit_processing.sh reco|ntuple TAG [--test]
set -euo pipefail

usage() {
  echo "Usage: $0 premix TAG GENSIM_DIR MINBIAS_DIR [--test]" >&2
  echo "       $0 digiraw TAG GENSIM_DIR [--test]" >&2
  echo "       $0 reco|ntuple TAG [--test]" >&2
  exit 2
}

(( $# >= 2 )) || usage
stage=$1
tag=$2
[[ $tag =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || usage
shift 2
test_mode=0
if (( $# > 0 )) && [[ ${@: -1} == --test ]]; then
  test_mode=1
  remaining=$(( $# - 1 ))
  set -- "${@:1:remaining}"
fi

store="$HOME/workspace/.store/deepmuonreco/chunk/Run3/MC"
output_tag=$tag
if (( test_mode )); then output_tag="$tag-test"; fi
secondary_dir=
ratio=0
case $stage in
  premix)
    (( $# == 2 )) || usage
    source_stage=gensim
    input_dir=$1
    secondary_stage=minbias
    secondary_dir=$2
    ratio=10
    if (( test_mode )); then ratio=1; fi
    ;;
  digiraw)
    (( $# == 1 )) || usage
    source_stage=gensim
    input_dir=$1
    secondary_stage=premix
    secondary_dir=$store/$output_tag/premix
    ;;
  reco|ntuple)
    (( $# == 0 )) || usage
    if [[ $stage == reco ]]; then source_stage=digiraw; else source_stage=reco; fi
    input_dir=$store/$output_tag/$source_stage
    ;;
  *) usage ;;
esac

check_series() {
  local directory=$1 prefix=$2 i number expected
  [[ -d $directory ]] || { echo "Missing $prefix directory: $directory" >&2; exit 1; }
  series_dir=$(cd -- "$directory" && pwd -P)
  if (( test_mode )); then
    expected="$series_dir/${prefix}_00000.root"
    [[ -s $expected ]] || { echo "Missing or empty $prefix input: $expected" >&2; exit 1; }
    series_files=("$expected")
    return
  fi
  shopt -s nullglob
  series_files=("$series_dir"/"${prefix}_"*.root)
  (( ${#series_files[@]} > 0 )) || {
    echo "No ${prefix}_*.root files in $series_dir" >&2
    exit 1
  }
  for ((i=0; i<${#series_files[@]}; i++)); do
    printf -v number '%05d' "$i"
    expected="$series_dir/${prefix}_$number.root"
    [[ ${series_files[i]} == "$expected" && -s $expected ]] || {
      echo "Missing or empty $prefix input: $expected" >&2
      exit 1
    }
  done
}

check_series "$input_dir" "$source_stage"
input_dir=$series_dir
inputs=("${series_files[@]}")
if [[ -n $secondary_dir ]]; then
  check_series "$secondary_dir" "$secondary_stage"
  secondary_dir=$series_dir
  secondaries=("${series_files[@]}")
  if (( test_mode )); then
    required=1
  elif [[ $stage == premix ]]; then
    required=$(( ${#inputs[@]} * ratio ))
  else
    required=${#inputs[@]}
  fi
  (( ${#secondaries[@]} >= required )) || {
    echo "Incomplete $secondary_stage pool: need $required files, found ${#secondaries[@]}" >&2
    exit 1
  }
  echo "$secondary_stage inputs: ${#secondaries[@]} files"
fi

here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd -- "$here/../.." && pwd)
count=${#inputs[@]}
events=-1
if [[ $stage == premix ]]; then events=1000; fi
if (( test_mode )); then count=1; events=1; fi
output="$store/$output_tag/$stage"
logs="$store/$output_tag/logs/$stage"
echo "$stage: $count jobs from $input_dir"
echo "output: $output"
mkdir -p "$logs" "$output"

for ((offset=0; offset<count; offset+=10000)); do
  batch=$(( count - offset ))
  (( batch <= 10000 )) || batch=10000
  submitted=$(sbatch --parsable --job-name="run3-$stage" \
    --array="0-$((batch - 1))" --output="$logs/%A_%a.out" \
    "$here/run_processing_array.sbatch" \
    "$stage" "$events" "$output" "$offset" "$input_dir" "$secondary_dir" "$ratio" "$repo")
  echo "submitted tasks $offset-$((offset + batch - 1)): job ${submitted%%;*}"
done
