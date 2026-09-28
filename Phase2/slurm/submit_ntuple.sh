#!/usr/bin/env bash
# Usage: bash Phase2/slurm/submit_ntuple.sh TAG [--test]
set -euo pipefail

if (( $# < 1 || $# > 2 )) || { (( $# == 2 )) && [[ $2 != --test ]]; }; then
  echo "Usage: $0 TAG [--test]" >&2
  exit 2
fi

tag=$1
[[ $tag =~ ^[a-zA-Z0-9][a-zA-Z0-9._-]*$ ]] || {
  echo "TAG must contain only letters, digits, dots, underscores, and hyphens" >&2
  exit 2
}
[[ -f $HOME/workspace/deepmuonreco/CMSSW_14_0_9/.deepmuonreco-build-status ]] &&
  [[ $(< "$HOME/workspace/deepmuonreco/CMSSW_14_0_9/.deepmuonreco-build-status") == ready ]] || {
    echo "Build CMSSW first: bash Phase2/setup_cmssw.sh" >&2
    exit 1
  }

max_events=-1
if [[ ${2:-} == --test ]]; then
  tag="$tag-test"
  max_events=1
fi

here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
sample="$HOME/workspace/.store/deepmuonreco/chunk/Phase2/MC/$tag"
rereco="$sample/rereco"
[[ -r $rereco/input-files.txt ]] || {
  echo "No ReReco input list for $tag; run submit_rereco.sh first" >&2
  exit 1
}
mapfile -t source_inputs < "$rereco/input-files.txt"
(( ${#source_inputs[@]} > 0 )) || { echo "ReReco input list is empty" >&2; exit 1; }

inputs=()
for ((i=0; i<${#source_inputs[@]}; i++)); do
  printf -v number '%05d' "$i"
  file="$rereco/rereco_$number.root"
  [[ -s $file ]] || { echo "Missing ReReco result: $file" >&2; exit 1; }
  inputs+=("$file")
done

output="$sample/ntuple"
logs="$sample/logs/ntuple"
mkdir -p "$output" "$logs"
manifest="$output/input-files.txt"
if [[ -e $manifest ]]; then
  cmp -s "$manifest" <(printf '%s\n' "${inputs[@]}") || {
    echo "The ReReco inputs for $tag changed; choose a new TAG" >&2
    exit 1
  }
else
  printf '%s\n' "${inputs[@]}" > "$manifest"
fi

echo "Ntuple: ${#inputs[@]} ReReco files -> $output"
for ((offset=0; offset<${#inputs[@]}; offset+=10000)); do
  count=$(( ${#inputs[@]} - offset ))
  (( count <= 10000 )) || count=10000
  jobid=$(sbatch --parsable --job-name=phase2-ntuple \
    --array="0-$((count - 1))" --output="$logs/%A_%a.out" \
    "$here/run_phase2.sbatch" \
    ntuple "$manifest" "$output" "$max_events" "$offset")
  jobid=${jobid%%;*}
  echo "Submitted inputs $offset-$((offset + count - 1)): $jobid"
done
