#!/usr/bin/env python3
import re
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

import util
from util import BENCHMARK_DIRECTORY


def extract_timestamp(benchmark_dir, threads=10, seqs=1000000):
    df = pd.read_csv(benchmark_dir / "performance_overview.csv")
    if threads:
        df = df[df["threads"] == threads]
    if seqs:
        df = df[df["seqs"] == seqs]
    df["benchmark"] = benchmark_dir.name
    return df


def main():
    timestamps = pd.concat(extract_timestamp(bm_dir, seqs=None)
                  for bm_dir in Path().iterdir()
                  if bm_dir.name.startswith("benchmark_"))
    timestamps.pivot(index="benchmark", columns=["seqs"], values="nucleotide alignment").T.plot(ylabel="time of nucleotide ali. in s", xlabel="Number of sequences in sample")
    plt.show()
    print(timestamps)

if __name__ == '__main__':
    main()
