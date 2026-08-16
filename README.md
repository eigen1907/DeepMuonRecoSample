# DeepMuonReco sample

```text
Ntuplizer/  shared plugin and notebook
Phase2/     Phase-2 configs and Condor submission
Run3/       Run 3 configs and Condor submission
```

## Phase-2

Use `CMSSW_14_0_9`.

```bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /path/to/CMSSW_14_0_9/src
cmsenv
scram b -j 8

cmsRun "$CMSSW_BASE/src/DeepMuonRecoSample/Phase2/test/runReReco_cfg.py" \
  inputFiles=file:/path/to/input.root outputFile=rereco.root maxEvents=-1

cmsRun "$CMSSW_BASE/src/DeepMuonRecoSample/Phase2/test/runDeepMuonRecoNtuplizer_cfg.py" \
  inputFiles=file:rereco.root outputFile=ntuple.root maxEvents=-1
```

The Phase-2 physics configuration and ntuple schema are unchanged; only the
paths moved under `Phase2/`.

### Phase-2 Condor

The same submit file is used for ReReco and ntuple jobs. Each line of the input
list is one absolute ROOT file path.

```bash
cd "$CMSSW_BASE/src/DeepMuonRecoSample/Phase2/condor"
mkdir -p logs /hdfs/path/to/rereco
find /hdfs/path/to/official-input -type f -name '*.root' -size +1000c \
  | sort > inputs-rereco.txt

condor_submit -batch-name phase2-rereco \
  cfg=runReReco_cfg.py label=rereco max_events=-1 \
  input_list="$PWD/inputs-rereco.txt" output_dir=/hdfs/path/to/rereco \
  submit.sub
```

Submit and worker nodes must see the same CMSSW project path.

After ReReco finishes, make a new input list and change only the cfg and output
directory:

```bash
mkdir -p /hdfs/path/to/ntuples
find /hdfs/path/to/rereco -type f -name '*.root' | sort > inputs-ntuple.txt

condor_submit -batch-name phase2-ntuple \
  cfg=runDeepMuonRecoNtuplizer_cfg.py label=ntuple max_events=-1 \
  input_list="$PWD/inputs-ntuple.txt" output_dir=/hdfs/path/to/ntuples \
  submit.sub
```

## Run 3

Use `CMSSW_14_0_21_patch1` with the package linked into its `src` directory.

```bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
scram project CMSSW CMSSW_14_0_21_patch1
cd CMSSW_14_0_21_patch1/src
ln -s /absolute/path/to/DeepMuonRecoSample DeepMuonRecoSample
cmsenv
scram b -j 8
```

Change to an empty working directory, then run the five configs:

```bash
mkdir -p /path/to/run3-work
cd /path/to/run3-work

cmsRun "$CMSSW_BASE/src/DeepMuonRecoSample/Run3/test/runGENSIM_cfg.py"
cmsRun "$CMSSW_BASE/src/DeepMuonRecoSample/Run3/test/runMinBiasGENSIM_cfg.py"
cmsRun "$CMSSW_BASE/src/DeepMuonRecoSample/Run3/test/runDIGIRAW_cfg.py"
cmsRun "$CMSSW_BASE/src/DeepMuonRecoSample/Run3/test/runRECO_cfg.py"
cmsRun "$CMSSW_BASE/src/DeepMuonRecoSample/Run3/test/runDeepMuonRecoNtuplizer_cfg.py"
```

The default files are `gensim.root`, `minbias.root`, `digiraw.root`,
`reco.root`, and `ntuple.root`. The Run 3 setup uses the official
`SingleMuFlatPt2To100` gun and generates TuneCP5 13.6 TeV MinBias events with
the local fragment. DIGI then reads those MinBias files for direct Summer24
pileup mixing.

Run 3 deliberately uses the same ntuple schema as Phase-2. Consequently,
`track_is_matched_muon` does not distinguish a signal muon from a pileup muon.

The default 1000-event MinBias sample is only for a local smoke test. For
another job, use a separate directory and pass the same `jobIndex=N` to the
first three configs.

### Run 3 Condor

Run 3 is split into separate jobs and keeps every stage output:

```text
signal GENSIM ─┐
               ├─> DIGIRAW ─> RECO ─> Ntuple
MinBias pool ──┘
```

Submit each stage only after checking that the previous stage finished and its
ROOT files were transferred successfully. Start from a built
`CMSSW_14_0_21_patch1` environment:

```bash
cd "$CMSSW_BASE/src/DeepMuonRecoSample/Run3/condor"
campaign_dir=/hdfs/path/to/run3
mkdir -p logs "$campaign_dir"/{gensim,minbias,digiraw,reco,ntuple}
```

Generate the signal and the shared MinBias pool independently. The ranges are
inclusive; choose the MinBias range and events per job for the campaign size.

```bash
condor_submit -batch-name run3-gensim \
  stage=gensim first_job=0 last_job=99 events=10 \
  output_dir="$campaign_dir/gensim" generation.sub

condor_submit -batch-name run3-minbias \
  stage=minbias first_job=0 last_job=99 events=1000 \
  output_dir="$campaign_dir/minbias" generation.sub
```

After both productions finish, create the shared pileup list and the indexed
signal input list. The pileup list contains CMSSW file names; the indexed list
contains `JOB_INDEX ABSOLUTE_PATH` on each line.

```bash
find "$campaign_dir/minbias" -type f -name 'minbias_*.root' -size +1000c \
  | sort | sed 's#^#file:#' > inputs-minbias.txt

find "$campaign_dir/gensim" -type f -name 'gensim_*.root' -size +1000c \
  | sed -E 's#.*/gensim_([0-9]+)\.root#\1 &#' \
  | sort -n > inputs-gensim.txt

condor_submit -batch-name run3-digiraw \
  stage=digiraw events=10 input_list="$PWD/inputs-gensim.txt" \
  pileup_list="$PWD/inputs-minbias.txt" \
  output_dir="$campaign_dir/digiraw" processing.sub
```

Build the next indexed list after each stage succeeds:

```bash
find "$campaign_dir/digiraw" -type f -name 'digiraw_*.root' -size +1000c \
  | sed -E 's#.*/digiraw_([0-9]+)\.root#\1 &#' \
  | sort -n > inputs-digiraw.txt

condor_submit -batch-name run3-reco \
  stage=reco events=-1 input_list="$PWD/inputs-digiraw.txt" \
  output_dir="$campaign_dir/reco" processing.sub

find "$campaign_dir/reco" -type f -name 'reco_*.root' -size +1000c \
  | sed -E 's#.*/reco_([0-9]+)\.root#\1 &#' \
  | sort -n > inputs-reco.txt

condor_submit -batch-name run3-ntuple \
  stage=ntuple events=-1 input_list="$PWD/inputs-reco.txt" \
  output_dir="$campaign_dir/ntuple" processing.sub
```

The retained files are `gensim/gensim_INDEX.root`,
`minbias/minbias_INDEX.root`, `digiraw/digiraw_INDEX.root`,
`reco/reco_INDEX.root`, and `ntuple/ntuple_INDEX.root`. Do not assign one
signal job index to two different event chunks; retrying the same failed job
should reuse its index. Submit and worker nodes must see the same CMSSW project
and mounted HDFS paths.
