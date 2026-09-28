#!/usr/bin/env bash
# Usage: bash submit_gensim.sh [--test|--minbias-test|--minbias JOBS]
set -euo pipefail

case ${1:-} in
  --test) stage=gensim; events=1; jobs=1; tag=run3-gensim-smoke-v002; time=(--time=02:00:00) ;;
  '') stage=gensim; events=1000; jobs=249; tag=run3-singlemu-gensim-v002; time=() ;;
  --minbias-test) stage=minbias; events=100; jobs=1; tag=run3-minbias-smoke-v001; time=(--time=04:00:00) ;;
  --minbias) stage=minbias; events=1000; jobs=${2:-}; tag=run3-minbias-v001; time=() ;;
  *) echo "Usage: $0 [--test|--minbias-test|--minbias JOBS]" >&2; exit 2 ;;
esac
if [[ ! $jobs =~ ^[1-9][0-9]*$ ]] || \
   { [[ ${1:-} == --minbias ]] && (( $# != 2 )); } || \
   { [[ ${1:-} != --minbias ]] && (( $# > 1 )); }; then
  echo "Invalid arguments or job count." >&2
  exit 2
fi

here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
store="$HOME/workspace/.store/deepmuonreco/chunk/Run3/MC"
output="$store/$tag/$stage"
logs="$store/$tag/logs/$stage"
mem=8G
if [[ $stage == minbias ]]; then mem=3G; fi
echo "$jobs jobs × $events events = $((jobs * events)) events"
echo "output: $output"
mkdir -p "$logs" "$output"
for ((offset=0; offset<jobs; offset+=10000)); do
  count=$(( jobs - offset ))
  (( count <= 10000 )) || count=10000
  submitted=$(sbatch --parsable --job-name="run3-$stage" \
    --array="0-$((count - 1))" --mem="$mem" \
    --output="$logs/%A_%a.out" \
    "${time[@]}" "$here/run_gensim_array.sbatch" "$stage" "$events" "$output" "$offset")
  echo "submitted tasks $offset-$((offset + count - 1)): job ${submitted%%;*}"
done
