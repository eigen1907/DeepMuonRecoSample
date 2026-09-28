#!/usr/bin/env bash
# Usage: bash Phase2/slurm/submit_rereco.sh INPUT_LIST TAG [--test]
set -euo pipefail

if (( $# < 2 || $# > 3 )) || { (( $# == 3 )) && [[ $3 != --test ]]; }; then
  echo "Usage: $0 INPUT_LIST TAG [--test]" >&2
  exit 2
fi

input_list=$1
tag=$2
[[ $tag =~ ^[a-zA-Z0-9][a-zA-Z0-9._-]*$ ]] || {
  echo "TAG must contain only letters, digits, dots, underscores, and hyphens" >&2
  exit 2
}
[[ -r $input_list ]] || { echo "Cannot read $input_list" >&2; exit 1; }
[[ -f $HOME/workspace/deepmuonreco/CMSSW_14_0_9/.deepmuonreco-build-status ]] &&
  [[ $(< "$HOME/workspace/deepmuonreco/CMSSW_14_0_9/.deepmuonreco-build-status") == ready ]] || {
    echo "Build CMSSW first: bash Phase2/setup_cmssw.sh" >&2
    exit 1
  }

mapfile -t inputs < <(awk 'NF && $1 !~ /^#/ { print $1 }' "$input_list")
(( ${#inputs[@]} > 0 )) || { echo "Input list is empty" >&2; exit 2; }
for input in "${inputs[@]}"; do
  case "$input" in
    root://*|/store/*|/*|file:/*) ;;
    *) echo "Expected an absolute path, /store path, or root:// URL: $input" >&2; exit 2 ;;
  esac
done

max_events=-1
if [[ ${3:-} == --test ]]; then
  inputs=("${inputs[0]}")
  tag="$tag-test"
  max_events=1
fi

here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
store="$HOME/workspace/.store/deepmuonreco/chunk/Phase2/MC"
output="$store/$tag/rereco"
logs="$store/$tag/logs/rereco"
mkdir -p "$output" "$logs"
manifest="$output/input-files.txt"
if [[ -e $manifest ]]; then
  cmp -s "$manifest" <(printf '%s\n' "${inputs[@]}") || {
    echo "The input list for $tag changed; choose a new TAG" >&2
    exit 1
  }
else
  printf '%s\n' "${inputs[@]}" > "$manifest"
fi

echo "ReReco: ${#inputs[@]} files -> $output"
for ((offset=0; offset<${#inputs[@]}; offset+=10000)); do
  count=$(( ${#inputs[@]} - offset ))
  (( count <= 10000 )) || count=10000
  jobid=$(sbatch --parsable --job-name=phase2-rereco \
    --array="0-$((count - 1))" --output="$logs/%A_%a.out" \
    "$here/run_phase2.sbatch" \
    rereco "$manifest" "$output" "$max_events" "$offset")
  jobid=${jobid%%;*}
  echo "Submitted inputs $offset-$((offset + count - 1)): $jobid"
done
