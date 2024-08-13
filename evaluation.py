#!/usr/bin/env python3
import re

import pandas as pd

import util
from util import BENCHMARK_DIRECTORY


def extract_timestamps(benchmark_directory):
    print("Reading logs from", benchmark_directory.name)
    threads, seqs = re.fullmatch(r"benchmark_(\d+)threads_(\d+)seqs", benchmark_directory.name).groups()
    log = util.read_benchmark_log_file(threads, seqs)

    timestamps = re.finditer("INFO: TIMESTAMP: Completed\s+([A-Za-z ]+)\s+:\s+(\d+)\s+seconds", log)
    data = {name.strip(): int(seconds) for name, seconds in (ts.groups() for ts in timestamps)}
    data["threads"] = int(threads)
    data["seqs"] = int(seqs)
    return data

def main():
    timestamps = (extract_timestamps(bm_dir) for bm_dir in BENCHMARK_DIRECTORY.iterdir())
    timestamps = pd.DataFrame.from_records(timestamps)
    timestamps.set_index(["seqs", "threads"], inplace=True, drop=True)
    timestamps.sort_index(inplace=True)
    timestamps.to_csv(str(BENCHMARK_DIRECTORY / "performance_overview.csv"))
    print(timestamps)

if __name__ == '__main__':
    main()