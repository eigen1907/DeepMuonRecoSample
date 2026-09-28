# Run 3 Slurm 작업

`DeepMuonRecoSample` 저장소에서 명령을 실행한다. 한 단계를 제출하고 결과 파일이 모두 완성된 뒤 다음 단계를 제출한다.

`signal GEN-SIM + MinBias GEN-SIM → premix → DIGI-RAW → RECO → MC ntuple`

새 MC 경로는 CMSSW premix를 사용한다. 기존 `Run3/test/run_mc_digiraw_cfg.py`와 Condor 작업은 MinBias를 DIGI-RAW에서 직접 섞으므로 별도 경로로 둔다.

## CMSSW 한 번만 준비

GEN-SIM, premix, DIGI-RAW, RECO는 EL8 Apptainer 안에서 CVMFS의 `CMSSW_14_0_21_patch1`을 사용한다. ntuple 플러그인은 `~/workspace/deepmuonreco/CMSSW_14_0_21_patch1`에 한 번 빌드해야 한다. Apptainer가 있는 lugia 셸에서 실행한다.

```bash
cd ~/workspace/deepmuonreco/DeepMuonRecoSample
bash Run3/setup_cmssw.sh
```

빌드에 성공하면 `DONE`이 출력된다. ntuple Slurm 작업은 빌드 상태가 `ready`인지 확인한다.

## 현재 MC 생산에 이어서 작업하기

이미 제출한 signal 및 MinBias 배열은 **기존의 평평한 디렉터리**에 계속 출력한다. premix 입력으로 그대로 사용하면 되며, 파일을 옮기거나 다시 제출할 필요가 없다.

```bash
SIG="$HOME/workspace/deepmuonreco/run3-slurm-output/run3-singlemu-gensim-v002"
MB="$HOME/workspace/deepmuonreco/run3-slurm-output/run3-minbias-v001"
TAG=run3-singlemu-premix-v001
```

제출된 작업이 사용하는 경로는 `scontrol show job <job-id>`로도 확인할 수 있다. **새로** 제출하는 작업의 출력 경로는 `/users/hep/joshin/workspace/.store/deepmuonreco/chunk/Run3/MC/<tag>/<stage>/`이다. 따라서 `submit_gensim.sh --minbias 2490`을 다시 실행하면 새 저장 경로에 MinBias 생산을 하나 더 시작하게 된다.

현재 signal은 249개 파일(`gensim_00000.root`부터 `gensim_00248.root`), MinBias 목표는 2490개 파일이다. 생산용 premix 작업 `i`는 `minbias_(10i)`부터 `minbias_(10i+9)`까지 **MinBias 파일 10개**를 읽고 premix 이벤트 1000개를 만든다. 제출 스크립트는 MinBias 파일이 1:10 비율로 모두 준비됐는지 확인한다. DIGI-RAW 작업 `i`는 signal 파일 `i`와 premix 파일 `i`를 합친다.

## 1 이벤트 시험 후 전체 제출

`gensim_00000.root`와 `minbias_00000.root`가 준비되면 아래 명령으로 단계별 1 이벤트 시험을 한다. **각 단계가 끝나고 결과 파일을 확인한 다음** 다음 명령을 실행한다.

```bash
bash Run3/slurm/submit_processing.sh premix "$TAG" "$SIG" "$MB" --test
bash Run3/slurm/submit_processing.sh digiraw "$TAG" "$SIG" --test
bash Run3/slurm/submit_processing.sh reco "$TAG" --test
bash Run3/slurm/submit_processing.sh ntuple "$TAG" --test
```

시험 결과는 `~/workspace/.store/deepmuonreco/chunk/Run3/MC/${TAG}-test/<stage>/`에 생성되므로 생산 결과와 섞이지 않는다. 각 단계가 끝난 뒤 로그에 `DONE`이 있는지 확인한다. 새 MC 로그 경로는 같은 `<tag>/logs/<stage>/`이다.

**MinBias 2490개가 모두 완성되고** 시험 체인이 성공하면 아래 순서로 전체 생산을 제출한다. 단계마다 전체 결과 파일이 나온 뒤 다음 명령을 실행한다.

```bash
bash Run3/slurm/submit_processing.sh premix "$TAG" "$SIG" "$MB"
bash Run3/slurm/submit_processing.sh digiraw "$TAG" "$SIG"
bash Run3/slurm/submit_processing.sh reco "$TAG"
bash Run3/slurm/submit_processing.sh ntuple "$TAG"
```

배열의 동시 실행 수는 Slurm의 가용 자원에 따른다. 대기 및 실행 상태는 `squeue -u "$USER"`로 확인한다. 결과 파일명은 `<stage>_00000.root`, `<stage>_00001.root` 순서다. 같은 태그로 다시 제출할 때 내용이 있는 최종 결과 파일은 건너뛴다. 전체 MC 결과는 `~/workspace/.store/deepmuonreco/chunk/Run3/MC/$TAG/<stage>/`에 저장된다.

새 생산에 필요한 GEN-SIM 및 MinBias 제출 명령은 다음과 같다. 이 명령들 역시 **새 저장 경로**를 사용한다.

```bash
bash Run3/slurm/submit_gensim.sh --test
bash Run3/slurm/submit_gensim.sh --minbias-test
bash Run3/slurm/submit_gensim.sh
bash Run3/slurm/submit_gensim.sh --minbias 2490
```

현재 signal 및 MinBias 작업은 이미 제출되어 있으므로 위 두 생산 명령은 **새 생산을 시작할 때만** 실행한다. 새 배열에는 동시 실행 개수 제한이나 배열 간 의존성을 지정하지 않는다.
새로 생성한 GEN-SIM을 후속 단계의 입력으로 쓸 때는 다음 경로를 사용한다.

```bash
SIG="$HOME/workspace/.store/deepmuonreco/chunk/Run3/MC/run3-singlemu-gensim-v002/gensim"
MB="$HOME/workspace/.store/deepmuonreco/chunk/Run3/MC/run3-minbias-v001/minbias"
```

## Run 3 충돌 데이터 ntuple

충돌 데이터는 ntuple 단계만 필요하다. 같은 CMSSW 빌드를 마친 뒤 원격 입력을 읽을 수 있는 VOMS 프록시(`$X509_USER_PROXY` 또는 `~/.globus/cms-proxy`)를 준비한다. 기본 입력 목록은 선택된 Muon0 AOD 파일을 사용한다.

```bash
voms-setup
bash Run3/slurm/submit_data_ntuple.sh --test
bash Run3/slurm/submit_data_ntuple.sh
```

Data 결과는 `/users/hep/joshin/workspace/.store/deepmuonreco/chunk/Run3/Data/<production-id>/ntuple/`에 저장된다. Data 제출 스크립트는 `jobs.tsv`를 이 `ntuple` 디렉터리의 상위에, 로그를 같은 위치의 `logs/ntuple/`에 둔다.
