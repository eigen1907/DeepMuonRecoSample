# Phase2 MC: ReReco → ntuple

CMSSW 작업 영역은 **한 번만 직접** 준비합니다. Apptainer가 있는 `lugia`
셸에서 저장소 루트로 이동해 실행하세요.

```bash
bash Phase2/setup_cmssw.sh
```

플러그인 소스를 수정했다면 같은 명령으로 다시 빌드하세요.

## 1. MC DIGI-RAW → ReReco

입력 ROOT 파일 경로를 한 줄에 하나씩 담은 파일을 준비합니다. `/store/...root`,
`root://...`, 로컬 절대 경로를 사용할 수 있습니다. 원격 입력이면 작업 노드가
파일을 읽을 수 있도록 `~/.globus/cms-proxy`를 갱신해 두세요.

```bash
bash Phase2/slurm/submit_rereco.sh /path/to/digiraw-files.txt phase2-mc-v001 --test
bash Phase2/slurm/submit_rereco.sh /path/to/digiraw-files.txt phase2-mc-v001
```

`--test`는 첫 입력 파일에서 1 event만 처리하고 별도의
`phase2-mc-v001-test` 디렉터리에 저장합니다. 두 번째 명령은 목록의 모든
입력 파일을 처리합니다.

## 2. ReReco → ntuple

ReReco 작업이 모두 끝나면 같은 태그로 ntuple을 제출합니다. 스크립트가
ReReco 결과를 자동으로 찾아 입력 목록을 만들며, 빠진 결과가 있으면 제출을
멈춥니다.

```bash
bash Phase2/slurm/submit_ntuple.sh phase2-mc-v001 --test
bash Phase2/slurm/submit_ntuple.sh phase2-mc-v001
```

`ntuple --test`는 위의 `rereco --test` 결과가 생성된 후 실행합니다.

출력은 `/users/hep/joshin/workspace/.store/deepmuonreco/chunk/Phase2/MC/<tag>/rereco/` 및
`/users/hep/joshin/workspace/.store/deepmuonreco/chunk/Phase2/MC/<tag>/ntuple/`에 저장됩니다.
각 입력 파일당 Slurm 배열 작업 하나를 생성합니다. 동시 실행 수는 Slurm의 가용 자원에 따릅니다.
로그는 같은 `<tag>/logs/<stage>/`에 저장됩니다.

```bash
squeue -u "$USER"
tail -f "$HOME/workspace/.store/deepmuonreco/chunk/Phase2/MC/phase2-mc-v001/logs/rereco/JOBID_0.out"
```
