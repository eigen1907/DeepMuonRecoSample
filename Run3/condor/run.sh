#!/bin/bash
set -euo pipefail

stage="$1"
job_index="$2"
events="$3"
input_file="$4"
pileup_list="$5"

: "${DMR_CMSSW_BASE:?Missing DMR_CMSSW_BASE}"
: "${_CONDOR_SCRATCH_DIR:?Not running in a Condor scratch directory}"

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "${DMR_CMSSW_BASE}/src"
eval "$(scram runtime -sh)"

if [[ "${CMSSW_VERSION}" != "CMSSW_14_0_21_patch1" ]]; then
    echo "Run 3 requires CMSSW_14_0_21_patch1, got ${CMSSW_VERSION}" >&2
    exit 2
fi

cd "${_CONDOR_SCRATCH_DIR}"
cfg_dir="${DMR_CMSSW_BASE}/src/DeepMuonRecoSample/Run3/test"

case "${stage}" in
    gensim)
        cmsRun "${cfg_dir}/runGENSIM_cfg.py" \
            maxEvents="${events}" jobIndex="${job_index}" outputFile=output.root
        ;;
    minbias)
        cmsRun "${cfg_dir}/runMinBiasGENSIM_cfg.py" \
            maxEvents="${events}" jobIndex="${job_index}" outputFile=output.root
        ;;
    digiraw)
        pileup_list="${pileup_list##*/}"
        cmsRun "${cfg_dir}/runDIGIRAW_cfg.py" \
            inputFiles="file:${input_file}" secondaryInputList="${pileup_list}" \
            maxEvents="${events}" jobIndex="${job_index}" outputFile=output.root
        ;;
    reco)
        cmsRun "${cfg_dir}/runRECO_cfg.py" \
            inputFiles="file:${input_file}" maxEvents="${events}" outputFile=output.root
        ;;
    ntuple)
        cmsRun "${cfg_dir}/runDeepMuonRecoNtuplizer_cfg.py" \
            inputFiles="file:${input_file}" maxEvents="${events}" outputFile=output.root
        ;;
    *)
        echo "Unknown stage: ${stage}" >&2
        exit 2
        ;;
esac
