#!/usr/bin/env bash
# Usage: bash submit_gensim.sh [--test|--dry-run|--minbias-test|--minbias JOBS]
set -euo pipefail

case ${1:-} in
  --test) stage=gensim; events=1; jobs=1; parallel=1; tag=run3-gensim-smoke-v002; time=(--time=02:00:00) ;;
  ''|--dry-run) stage=gensim; events=1000; jobs=249; parallel=4; tag=run3-singlemu-gensim-v002; time=() ;;
  --minbias-test) stage=minbias; events=100; jobs=1; parallel=1; tag=run3-minbias-smoke-v001; time=(--time=04:00:00) ;;
  --minbias) stage=minbias; events=1000; jobs=${2:-}; parallel=4; tag=run3-minbias-v001; time=() ;;
  *) echo "Usage: $0 [--test|--dry-run|--minbias-test|--minbias JOBS]" >&2; exit 2 ;;
esac
if [[ ! $jobs =~ ^[1-9][0-9]*$ ]] || (( jobs > 100000 )) || \
   { [[ ${1:-} == --minbias ]] && (( $# != 2 )); } || \
   { [[ ${1:-} != --minbias ]] && (( $# > 1 )); }; then
  echo "Invalid arguments or job count (1–100000 jobs)." >&2
  exit 2
fi

here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
output="$HOME/workspace/deepmuonreco/run3-slurm-output/$tag"
logs="$here/logs/$tag"
echo "$jobs jobs × $events events = $((jobs * events)) events"
echo "output: $output"
dependency=()
for ((offset=0; offset<jobs; offset+=1000)); do
  count=$(( jobs - offset ))
  (( count > 1000 )) && count=1000
  cmd=(sbatch --parsable --job-name="run3-$stage"
    --array="0-$((count - 1))%$parallel" --output="$logs/%A_%a.out"
    "${time[@]}" "${dependency[@]}"
    "$here/run_gensim_array.sbatch" "$stage" "$events" "$output" "$offset")
  if [[ ${1:-} == --dry-run ]]; then
    printf 'command:'; printf ' %q' "${cmd[@]}"; printf '\n'
    continue
  fi
  mkdir -p "$logs" "$output"
  submitted=$("${cmd[@]}")
  jobid=${submitted%%;*}
  echo "submitted tasks $offset-$((offset + count - 1)): job $jobid"
  dependency=(--dependency="afterany:$jobid")
done
