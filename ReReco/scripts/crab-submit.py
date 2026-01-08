#!/usr/bin/env python3
"""
adapted from https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3AdvancedTutorial#Exercise_6_CRABAPI_library

source /cvmfs/cms.cern.ch/common/crab-setup.sh
"""
from http.client import HTTPException
import warnings
from pathlib import Path
import json
from typing import Optional
import getpass
from datetime import datetime
import argparse
from CRABClient.ClientExceptions import ClientException # type: ignore
from CRABClient.UserUtilities import config as CrabConfig # type: ignore
from CRABAPI.RawCommand import crabCommand # type: ignore


def submit(config: CrabConfig,
           input_dataset: str,
           lumi_mask: Optional[str],
           user: str,
           name: str
):
    # /PrimaryDataset/ProcessedDataset/DataTier
    _, primary, processed, _ = input_dataset.split('/')
    pset_stem = Path(config.JobType.psetName).stem

    config.Data.inputDataset = input_dataset
    if lumi_mask is not None:
        config.Data.lumiMask = lumi_mask
    config.Data.outLFNDirBase = f'/store/user/{user}/{name}'
    config.Data.outputDatasetTag = f'{name}'
    config.General.requestName = f'{primary}_{name}'[:100]

    try:
        print(f"Submitting for {input_dataset=}")
        crabCommand('submit', config=config)
    except HTTPException as error:
        warnings.warn(f"Submission for {input_dataset=} failed: {error}")
    except ClientException as error:
        warnings.warn(f"Submission for {input_dataset=} failed: {error}")

    return config


def run(pset: Path,
        input_list: list[dict[str, str]],
        user: str,
        storage_site: str,
        name: str,
):
    """
    """
    config = CrabConfig()
    # General
    config.General.transferLogs = False
    config.General.transferOutputs = True
    # JobType
    config.JobType.pluginName = 'Analysis'
    config.JobType.psetName = str(pset)
    # Data
    config.Data.publication  = False
    config.Data.allowNonValidInputDataset = True
    config.Data.splitting = 'FileBased'
    config.Data.unitsPerJob = 1
    #config.Data.splitting = 'EventAwareLumiBased'
    #config.Data.unitsPerJob = 12
    #config.Data.splitting = 'Automatic'
    # Site
    config.Site.storageSite = storage_site
    #config.Site.blacklist = ['T1_DE_KIT']
    # Set max cpus, memory (max 3 GB memory by 1 core by default)
    config.JobType.numCores = 8
    config.JobType.maxMemoryMB = 16000

    for item in input_list:
        submit(config=config, input_dataset=item['input_dataset'],
               lumi_mask=item.get('lumi_mask', None), user=user, name=name)


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('-p', '--pset', type=Path, help='cfg file')
    parser.add_argument('-i', '--input', type=Path,
                        help='json file containing input datasets')
    parser.add_argument('-s', '--storage-site', default='T3_CH_CERNBOX',
                        type=str, help='storage site name')
    parser.add_argument('-u', '--user', default=getpass.getuser(), type=str,
                        help='lxplus user id')
    parser.add_argument('-n', '--name', default='RPC-TnP-NanoAOD', type=str,
                        help='project name')
    args = parser.parse_args()

    if not args.pset.exists():
        raise FileNotFoundError(f'{args.pset}')

    if not args.input.exists():
        raise FileNotFoundError(f'{args.input}')

    with open(args.input) as stream:
        input_list = json.load(stream)

    run(
        pset=args.pset,
        input_list=input_list,
        user=args.user,
        storage_site=args.storage_site,
        name=args.name
    )


if __name__ == '__main__':
    main()