#!/usr/bin/env python3
"""Split listed data files into 1,000-event jobs and submit them."""

import argparse
import getpass
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys


CHUNK_SIZE = 1000
MAX_EVENTS = 250_000
PRODUCTION_ID = re.compile(r"^[a-z0-9][a-z0-9.-]*$")


def die(message):
    print(f"error: {message}", file=sys.stderr)
    raise SystemExit(2)


def run(command):
    return subprocess.run(
        command,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )


def resolve_file_list(selector, repo_root):
    path = Path(selector).expanduser()
    if not path.exists():
        path = repo_root / "Run3" / "inputs" / f"{selector}.txt"
    if not path.is_file() or path.suffix != ".txt":
        die(f"file list not found: {selector}")
    if not PRODUCTION_ID.fullmatch(path.stem):
        die("file-list name must use lowercase letters, numbers, dots, and hyphens")
    return path.resolve(), path.stem


def read_inputs(path):
    inputs = []
    for line_number, raw_line in enumerate(path.read_text().splitlines(), 1):
        value = raw_line.split("#", 1)[0].strip()
        if not value:
            continue
        if any(character.isspace() for character in value):
            die(f"{path}:{line_number}: put one ROOT path on each line")
        if value.startswith("/store/"):
            value = "root://cms-xrd-global.cern.ch/" + value
        elif value.startswith("/"):
            value = "file:" + value
        elif not value.startswith(("root://", "file:")):
            die(f"{path}:{line_number}: invalid ROOT path")
        if not value.endswith(".root"):
            die(f"{path}:{line_number}: input must end with .root")
        if value in inputs:
            die(f"duplicate input: {value}")
        inputs.append(value)
    if not inputs:
        die(f"file list is empty: {path}")
    return inputs


def require_environment(repo_root):
    if (os.environ.get("CMSSW_VERSION") != "CMSSW_14_0_21_patch1"
            or not os.environ.get("CMSSW_BASE")):
        die("run cmsenv in CMSSW_14_0_21_patch1 first")
    cmssw_base = Path(os.environ["CMSSW_BASE"]).resolve()
    package = cmssw_base / "src" / "DeepMuonRecoSample"
    if package.resolve() != repo_root.resolve():
        die(f"CMSSW uses a different package: {package.resolve()}")

    proxy_result = run(["voms-proxy-info", "-path"])
    if proxy_result.returncode != 0:
        die("run: voms-proxy-init --voms cms --valid 192:00")
    proxy = Path(proxy_result.stdout.strip())
    valid = subprocess.run(
        ["voms-proxy-info", "-file", str(proxy), "-exists", "-valid", "1:00"],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    if valid.returncode != 0:
        die("CMS proxy is valid for less than one hour")
    return cmssw_base, proxy


def count_events(inputs):
    try:
        import ROOT
    except ModuleNotFoundError:
        die("PyROOT is unavailable; run cmsenv first")

    ROOT.gROOT.SetBatch(True)
    counts = []
    for index, input_file in enumerate(inputs, 1):
        print(f"counting   : {index}/{len(inputs)} {input_file}", flush=True)
        root_file = ROOT.TFile.Open(input_file)
        if not root_file or root_file.IsZombie():
            die(f"cannot open input: {input_file}")
        events_tree = root_file.Get("Events")
        events = int(events_tree.GetEntriesFast()) if events_tree else 0
        root_file.Close()
        if events <= 0:
            die(f"input has no Events entries: {input_file}")
        counts.append((input_file, events))
    return counts


def make_jobs(input_counts):
    chunks_by_file = []
    for input_file, file_events in input_counts:
        chunks_by_file.append([
            {
                "events": min(CHUNK_SIZE, file_events - skip_events),
                "skip": skip_events,
                "input": input_file,
            }
            for skip_events in range(0, file_events, CHUNK_SIZE)
        ])

    jobs = []
    selected_events = 0
    for chunk_index in range(max(map(len, chunks_by_file))):
        for chunks in chunks_by_file:
            if chunk_index < len(chunks):
                job = chunks[chunk_index].copy()
                job["events"] = min(job["events"], MAX_EVENTS - selected_events)
                job["output"] = f"ntuple-{len(jobs):04d}.root"
                jobs.append(job)
                selected_events += job["events"]
                if selected_events == MAX_EVENTS:
                    return jobs
    return jobs


def choose_paths(site, repo_root, production_id):
    work_root = Path(os.environ.get("DMR_WORK_ROOT", repo_root.parent / "prod"))
    if os.environ.get("DMR_OUTPUT_ROOT"):
        output_root = Path(os.environ["DMR_OUTPUT_ROOT"])
    elif site == "cern":
        user = getpass.getuser()
        output_root = Path(f"/eos/user/{user[0]}/{user}/deepmuonreco")
    else:
        die("set DMR_OUTPUT_ROOT on UOS")
    if not work_root.is_absolute() or not output_root.is_absolute():
        die("DMR_WORK_ROOT and DMR_OUTPUT_ROOT must be absolute paths")
    return work_root / production_id, output_root / production_id


def missing_jobs(output_dir, jobs):
    expected = {job["output"] for job in jobs}
    unexpected = [
        path for path in output_dir.glob("ntuple-*.root")
        if path.name not in expected
    ]
    if unexpected:
        die(f"unexpected output exists: {unexpected[0]}")

    complete = []
    missing = []
    for job in jobs:
        path = output_dir / job["output"]
        if not path.exists():
            missing.append(job)
        elif path.is_file() and path.stat().st_size >= 1024:
            complete.append(job)
        else:
            die(f"invalid output exists: {path}")
    return complete, missing


def active_jobs(production_id):
    result = run([
        "condor_q",
        "-constraint", f'JobBatchName == "{production_id}"',
        "-af", "ClusterId", "ProcId",
    ])
    if result.returncode != 0:
        die("cannot query Condor queue:\n" + result.stdout.strip())
    return [line for line in result.stdout.splitlines() if line.strip()]


def write_jobs(path, jobs):
    path.write_text("\n".join(
        f"{job['events']} {job['skip']} "
        f"{job['input']} {job['output']}"
        for job in jobs
    ) + "\n")


def submit(site, production_id, cmssw_base, proxy, work_dir, output_dir, jobs):
    queued = active_jobs(production_id)
    if queued:
        die(f"{len(queued)} jobs from this production are already in Condor")

    log_dir = work_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    persistent_proxy = work_dir / "x509up"
    shutil.copy2(proxy, persistent_proxy)
    persistent_proxy.chmod(0o600)
    job_list = work_dir / "jobs.tsv"
    write_jobs(job_list, jobs)

    if site == "cern":
        destination = f"root://eosuser.cern.ch/{output_dir}/"
        remaps = ""
    else:
        destination = ""
        remaps = f"$(output_file)={output_dir}/$(output_file)"

    condor_dir = Path(__file__).resolve().parent
    result = run([
        "condor_submit",
        "-batch-name", production_id,
        f"condor_dir={condor_dir}",
        f"cmssw_base={cmssw_base}",
        f"job_list={job_list}",
        f"log_root={log_dir}",
        f"proxy={persistent_proxy}",
        f"output_destination={destination}",
        f"output_remaps={remaps}",
        str(condor_dir / "submit_data.sub"),
    ])
    print(result.stdout, end="")
    if result.returncode != 0:
        raise SystemExit(result.returncode)


def main():
    parser = argparse.ArgumentParser(
        description="Count listed data events and submit 1,000-event jobs."
    )
    parser.add_argument(
        "production",
        nargs="?",
        default="data-run3-muon0-2024cde-v001",
        help="production ID (default: data-run3-muon0-2024cde-v001)",
    )
    parser.add_argument("--site", choices=("cern", "uos"), default="cern")
    parser.add_argument("--yes", action="store_true", help="submit without asking")
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[2]
    file_list, production_id = resolve_file_list(args.production, repo_root)
    inputs = read_inputs(file_list)
    cmssw_base, proxy = require_environment(repo_root)
    work_dir, output_dir = choose_paths(args.site, repo_root, production_id)
    input_counts = count_events(inputs)
    jobs = make_jobs(input_counts)
    complete, missing = missing_jobs(output_dir, jobs)

    print(f"production : {production_id}")
    print(f"input      : {len(inputs)} files, {sum(n for _, n in input_counts)} events")
    print(f"selected   : {sum(job['events'] for job in jobs)} events")
    print(f"jobs       : {len(jobs)} (<= 1,000 events/job)")
    print(f"output     : {output_dir}")
    print(f"status     : {len(complete)} complete, {len(missing)} missing")

    if not missing:
        print("nothing to submit")
        return
    if not args.yes:
        answer = input(f"Submit {len(missing)} jobs? [y/N] ").strip().lower()
        if answer not in {"y", "yes"}:
            print("not submitted")
            return
    submit(
        args.site, production_id, cmssw_base, proxy,
        work_dir, output_dir, missing,
    )


if __name__ == "__main__":
    main()
