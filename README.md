# DeepMuonRecoSample

DeepMuonReco 학습용 ROOT 샘플과 Ntuple을 만드는 CMSSW 패키지입니다.

```text
Ntuplizer/                     analyzer와 CMSSW Python 설정
Run3/
├── inputs/                    data ROOT 파일 목록
├── test/                      CMSSW 설정 파일
├── slurm/                     서버용 제출·실행 스크립트
└── condor/                    lxplus용 기존 작업
Phase2/
├── test/                      CMSSW 설정 파일
├── slurm/                     ReReco·ntuple 제출 스크립트
└── condor/                    기존 작업
```

`plugins/`, `python/`, `test/`는 CMSSW 표준 디렉터리 이름입니다.

## 서버에서 Slurm 실행

Apptainer가 있는 `lugia`에서 필요한 CMSSW 작업 영역을 **한 번만 직접**
준비합니다. Run3 GEN-SIM·premix·DIGI-RAW·RECO는 CVMFS의 CMSSW를 사용하며,
로컬 빌드는 Run3 ntuple·Data와 Phase2 작업에 필요합니다.

```bash
bash Phase2/setup_cmssw.sh   # Phase2 MC ReReco/ntuple용
bash Run3/setup_cmssw.sh     # Run3 MC/Data ntuple용
```

실제 Slurm 제출·실행 코드는 이 저장소의 `Run3/slurm/`, `Phase2/slurm/`에 있습니다.
단계별로 실행할 짧은 레시피는 형제 저장소
`../DeepMuonRecoStudies/runs/sample/`에 둡니다. 그 디렉터리의 `README.md`에
Run3 MC(`GEN-SIM + MinBias → premix → DIGI-RAW → RECO → ntuple`),
Run3 Data(ntuple), Phase2 MC(ReReco → ntuple)의 실행 순서가 있습니다.
각 단계의 입력 조건과 시험 방법은 [Run3](Run3/slurm/README.md),
[Phase2](Phase2/slurm/README.md) 문서를 참고하세요.

새 Slurm 결과와 로그는 다음 경로 아래에 저장됩니다. 이미 제출한 작업은
제출 당시의 출력 경로를 계속 사용합니다.

```text
~/workspace/.store/deepmuonreco/chunk/Run3/MC/<tag>/
~/workspace/.store/deepmuonreco/chunk/Run3/Data/<tag>/
~/workspace/.store/deepmuonreco/chunk/Phase2/MC/<tag>/
```

작업 상태는 `squeue -u "$USER"`, 배열 취소는 `scancel <job-id>`로 확인·처리합니다.

## 기존 lxplus Condor: Muon0 2024 CDE data

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

새로 생성하는 Data Ntuple에는 MC truth가 없으므로 `tp_*` 배열은 비어 있습니다.
트랙별 `track_is_matched_muon`, `track_match_tp_idx`, `track_match_quality`는
트랙 배열과 길이를 맞추되 모두 `-999`로 채웁니다. 재구성 정보인
`track_is_reco_muon`, `track_is_trk_muon`, `track_is_glb_muon`,
`track_is_pf_muon`은 Data에서도 실제 값을 기록합니다.

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

## UOS Condor

UOS에서는 worker가 볼 수 있는 경로만 지정하고 같은 명령을 사용합니다.

```bash
export DMR_WORK_ROOT=/path/to/deepmuonreco/prod
export DMR_OUTPUT_ROOT=/hdfs/your/path/deepmuonreco

./submit_data.py --site uos
```

## 기존 Condor MC와 Phase-2

기존 Condor Run3 MC는 `GENSIM + MinBias → DIGIRAW → RECO → ntuple`
순서로 직접 pileup mixing을 합니다. 관련 파일은
`Run3/condor/submit_mc_gensim.sub`, `submit_mc_processing.sub`,
`run_mc.sh`입니다. 새 Slurm 워크플로에는 별도 premix 단계가 있습니다.

Phase-2는 `CMSSW_14_0_9`에서 `Phase2/test/run_rereco_cfg.py` 실행 후
`Phase2/test/run_ntuple_cfg.py`를 실행합니다.
현재 서버에서 Apptainer와 Slurm으로 실행하는 방법은
[`Phase2/slurm/README.md`](Phase2/slurm/README.md)를 참고하세요.
