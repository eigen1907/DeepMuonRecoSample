#!/usr/bin/env python3
"""Print Run 3 Data ntuple tasks as events, skip, input, output TSV rows."""

import argparse
from pathlib import Path
import sys


CHUNK_SIZE = 1000
MAX_EVENTS = 250_000


def inputs_from(path):
    inputs = []
    for line_number, raw in enumerate(path.read_text().splitlines(), 1):
        value = raw.split("#", 1)[0].strip()
        if not value:
            continue
        if any(character.isspace() for character in value):
            raise ValueError(f"{path}:{line_number}: one ROOT path per line")
        if value.startswith("/store/"):
            value = "root://cms-xrd-global.cern.ch/" + value
        elif value.startswith("/"):
            value = "file:" + value
        elif not value.startswith(("root://", "file:")):
            raise ValueError(f"{path}:{line_number}: unsupported ROOT path")
        if not value.endswith(".root") or value in inputs:
            raise ValueError(f"{path}:{line_number}: invalid or duplicate ROOT path")
        inputs.append(value)
    if not inputs:
        raise ValueError(f"empty input list: {path}")
    return inputs


def event_count(input_file):
    import ROOT

    root_file = ROOT.TFile.Open(input_file)
    if not root_file or root_file.IsZombie():
        raise ValueError(f"cannot open input: {input_file}")
    try:
        tree = root_file.Get("Events")
        count = int(tree.GetEntriesFast()) if tree else 0
    finally:
        root_file.Close()
    if count <= 0:
        raise ValueError(f"input has no Events: {input_file}")
    return count


def jobs_for(input_counts):
    chunks_by_file = [
        [(min(CHUNK_SIZE, count - skip), skip, name)
         for skip in range(0, count, CHUNK_SIZE)]
        for name, count in input_counts
    ]
    selected = 0
    index = 0
    for chunk_number in range(max(map(len, chunks_by_file))):
        for chunks in chunks_by_file:
            if chunk_number >= len(chunks):
                continue
            events, skip, name = chunks[chunk_number]
            events = min(events, MAX_EVENTS - selected)
            yield events, skip, name, f"ntuple-{index:04d}.root"
            selected += events
            index += 1
            if selected == MAX_EVENTS:
                return


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_list", type=Path)
    parser.add_argument("--test", action="store_true",
                        help="plan one 10-event task without counting files")
    args = parser.parse_args()
    inputs = inputs_from(args.input_list)
    if args.test:
        print(f"10\t0\t{inputs[0]}\tntuple-0000.root")
        return

    input_counts = []
    for index, name in enumerate(inputs, 1):
        print(f"Counting {index}/{len(inputs)}: {name}", file=sys.stderr,
              flush=True)
        input_counts.append((name, event_count(name)))
    for events, skip, name, output in jobs_for(input_counts):
        print(f"{events}\t{skip}\t{name}\t{output}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as error:
        raise SystemExit(f"error: {error}") from error
