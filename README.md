# DeepMuonRecoSample

DeepMuonReco 학습용 Ntuple을 만드는 CMSSW 패키지입니다.

```text
Ntuplizer/                     analyzer와 확인 notebook
Run3/
├── inputs/                    data ROOT 파일 목록
├── test/                      CMSSW 설정 파일
└── condor/
    ├── submit_data.py         data 제출 명령
    ├── submit_data.sub        data Condor 설정
    ├── run_data.sh            data worker
    ├── submit_mc_gensim.sub   MC GENSIM 제출 설정
    ├── submit_mc_processing.sub
    └── run_mc.sh              MC worker
Phase2/                        Phase-2 설정과 Condor 코드
```

`plugins/`, `python/`, `test/`는 CMSSW 표준 디렉터리 이름입니다.

## Muon0 2024 CDE data

입력 목록은 ROOT 경로를 한 줄에 하나씩 적습니다.

```text
Run3/inputs/data-run3-muon0-2024cde-v001.txt
```

현재 목록에는 Run2024C/D/E AOD 파일 14개가 들어 있습니다. 다른 입력을
처리하려면 목록을 복사하고 `v002`처럼 production 버전을 올린 뒤
`/store/...root` 경로를 한 줄씩 적습니다. event 수는 작성하지 않습니다.

환경과 인증을 준비합니다.

```bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/user/j/joshin/workspace/deepmuonreco/CMSSW_14_0_21_patch1/src
cmsenv
scram b -j 8

kinit "$USER@CERN.CH"
voms-proxy-init --voms cms --valid 192:00
```

아래 명령 하나가 각 파일의 event 수를 읽고, 최대 250,000 events를 1,000
events씩 나눈 뒤 요약을 보여줍니다. 입력이 250,000 events보다 적으면 있는
만큼 전부 사용합니다. `y`를 입력해야 실제로 제출됩니다.

```bash
cd DeepMuonRecoSample/Run3/condor
./submit_data.py
```

출력 이름은 job 순서대로 `ntuple-0000.root`, `ntuple-0001.root`, ...가 됩니다.
입력 파일과 event 구간은 작업 디렉터리의 `jobs.tsv`에 기록됩니다.

작업과 출력 경로는 production ID로 연결됩니다.

```text
/afs/cern.ch/user/j/joshin/workspace/deepmuonreco/prod/<production-id>/
/eos/user/j/joshin/deepmuonreco/<production-id>/
```

같은 명령을 다시 실행하면 다음과 같이 동작합니다.

- Condor에 같은 production이 남아 있으면 제출하지 않습니다.
- queue가 끝났고 ROOT 파일이 빠졌다면 missing jobs만 다시 제출합니다.
- 모든 ROOT 파일이 있으면 아무것도 제출하지 않습니다.
- production의 입력 목록은 생성 후 변경하지 않습니다. 입력을 바꾸려면 production
  버전을 올려 새 목록을 만듭니다.

상태와 로그는 다음 위치에서 확인합니다.

```bash
condor_q -batch
condor_q -hold

ls /afs/cern.ch/user/j/joshin/workspace/deepmuonreco/prod/data-run3-muon0-2024cde-v001/logs
ls /eos/user/j/joshin/deepmuonreco/data-run3-muon0-2024cde-v001
```

기존 250,000-event 결과와 로그는 그대로 보존했습니다.

```text
/eos/user/j/joshin/deepmuonreco/Run3/Muon0_2024CDE/muon0_cde_v01/ntuple/
/afs/cern.ch/user/j/joshin/workspace/deepmuonreco/prod/run3-data-muon0-2024cde-250k-v001/logs/
```

## UOS

UOS에서는 worker가 볼 수 있는 경로만 지정하고 같은 명령을 사용합니다.

```bash
export DMR_WORK_ROOT=/path/to/deepmuonreco/prod
export DMR_OUTPUT_ROOT=/hdfs/your/path/deepmuonreco

./submit_data.py --site uos
```

## MC와 Phase-2

Run 3 MC는 `GENSIM + MinBias -> DIGIRAW -> RECO -> Ntuple` 순서이며
`Run3/condor/submit_mc_gensim.sub`, `submit_mc_processing.sub`,
`run_mc.sh`을 사용합니다.

Phase-2는 `CMSSW_14_0_9`에서 `Phase2/test/run_rereco_cfg.py` 실행 후
`Phase2/test/run_ntuple_cfg.py`를 실행합니다.

Ntuple 확인 notebook은 `Ntuplizer/notebooks/check-ntuple.ipynb`입니다.
