#!/usr/bin/env python
import sys
import os
import glob
from typing import Any, List, Dict
from pathlib import Path
from datetime import datetime
import argparse
import shutil
from htcondor.htcondor import Submit, Schedd # type: ignore
from coolname import generate_slug

def make_run_name(run_name: str | None) -> str:
    run_name = run_name or generate_slug(pattern=2)
    now = datetime.now().strftime('%y%m%d-%H%M%S')
    return f'{now}_{run_name}'

def get_root_files(input_dir: Path) -> List[Path]:
    if not input_dir.exists():
        raise FileNotFoundError(f"Directory not found: {input_dir}")
    
    files = list(input_dir.rglob('*.root'))
    
    valid_files = [f.absolute() for f in files if f.stat().st_size > 1000]
    return valid_files

def submit_jobs(
    pset: Path,
    input_dir: Path,
    file_list: List[Path],
    output_base: Path,
    run_name: str,
    cpus: int,
    memory: str,
    disk: str,
    max_events: int,
    node_list: List[str]
) -> None:
    dataset_label = input_dir.name
    job_batch_name = f'{run_name}.{dataset_label}'
    
    log_dir = Path('logs') / run_name / dataset_label
    log_dir.mkdir(parents=True, exist_ok=True)
    
    out_dir = output_base
    out_dir.mkdir(parents=True, exist_ok=True)

    executable = shutil.which('cmsRun')
    if not executable:
        raise FileNotFoundError('cmsRun not found. Did you run cmsenv?')

    pset_filename = pset.name
    pset_path = str(pset.absolute())

    requirements = '||'.join([f'(machine=="{each}.sscc.uos")' for each in node_list]) if node_list else 'True'
    
    job_items = []
    for idx, input_file in enumerate(file_list):
        output_filename = f"output_{idx}.root"
        
        job_items.append({
            'inputFile': f"file:{str(input_file)}",
            'finalOutput': str(out_dir / output_filename),
        })

    print(f" -> Found {len(job_items)} root files in {input_dir}")
    
    if not job_items:
        print("No files found. Exiting.")
        return

    submit_dict: Dict[str, Any] = {
        'universe': 'vanilla',
        'getenv': 'True',
        'executable': executable,
        
        'should_transfer_files': 'YES',
        'when_to_transfer_output': 'ON_EXIT',
        
        'transfer_input_files': pset_path,
        'transfer_output_files': 'temp.root',
        'transfer_output_remaps': '"temp.root = $(finalOutput)"',
        
        'arguments': f"{pset_filename} inputFiles=$(inputFile) outputFile=temp.root maxEvents={max_events}",

        'request_cpus': cpus,
        'request_memory': memory,
        'request_disk': disk,
        'JobBatchName': job_batch_name,
        'requirements': requirements,

        'output': str(log_dir / 'job_$(ClusterId).$(ProcId).out'),
        'error': str(log_dir / 'job_$(ClusterId).$(ProcId).err'),
        'log': str(log_dir / 'job_$(ClusterId).$(ProcId).log'),

        'stream_output': 'True',
        'stream_error': 'True',
    }

    sub = Submit(submit_dict)
    schedd = Schedd()
    
    submit_result = schedd.submit(sub, itemdata=iter(job_items))

    print(f'🚀 Submitted {len(job_items)} jobs. (Batch: {job_batch_name})')
    print(f'   Logs: {log_dir}')
    print(f'   Output: {out_dir}')


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Submit CMSSW jobs for a specific directory"
    )
    
    parser.add_argument('-p', '--pset', type=Path, required=True, help='CMSSW cfg file')
    
    parser.add_argument('-i', '--input', type=Path, required=True, dest='input_dir', help='Input directory containing root files')
    parser.add_argument('-o', '--output-dir', type=Path, default=Path.cwd() / 'results', help='Output directory')
    parser.add_argument('-n', '--maxEvents', type=int, default=-1, help='Number of events to process per job (-1 for all)')

    parser.add_argument('--run-name', type=str, default=None, help='Run name (optional)')
    parser.add_argument('--cpus', type=int, default=1, help='CPUs per job')
    parser.add_argument('--memory', type=str, default='8000MB', help='Memory per job')
    parser.add_argument('--disk', type=str, default='30GB', help='Disk space per job')

    parser.add_argument('--node', dest='node_list', type=str, nargs='*', default=[], help='Specific nodes')

    args = parser.parse_args()

    if not args.pset.exists():
        raise FileNotFoundError(f"PSet file not found: {args.pset}")
    if not args.input_dir.exists():
        raise FileNotFoundError(f"Input directory not found: {args.input_dir}")

    run_name = make_run_name(args.run_name)
    
    print("-" * 60)
    print(f"Input Dir : {args.input_dir}")
    print(f"Run Name  : {run_name}")
    print(f"Number of Events: {args.maxEvents}")
    print("-" * 60)

    files = get_root_files(args.input_dir)
    
    submit_jobs(
        pset=args.pset,
        input_dir=args.input_dir,
        file_list=files,
        output_base=args.output_dir,
        run_name=run_name,
        cpus=args.cpus,
        memory=args.memory,
        disk=args.disk,
        max_events=args.maxEvents,
        node_list=args.node_list
    )

if __name__ == "__main__":
    main()