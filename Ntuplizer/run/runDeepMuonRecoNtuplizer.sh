cmsenv
cmsRun ${CMSSW_BASE}/src/DeepMuonRecoSample/Ntuplizer/test/runDeepMuonRecoNtuplizer_cfg.py \
    inputFiles=file:input.root \
    outputFile=file:output.root \
    maxEvents=10