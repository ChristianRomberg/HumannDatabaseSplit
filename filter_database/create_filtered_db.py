#!/usr/bin/env python3
import bisect

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
from matplotlib.ticker import StrMethodFormatter

if __name__ == '__main__':
    cutoff=300000
    aa_counts = pd.read_csv("total_aa_counts.tsv", sep="\t", index_col=0, header=None, names=["protein", "count"])

    small_db_whitelist = set(aa_counts.index[:cutoff])

    uniref_input_db = "/home/chris/Projects/benchmark_HUMAnN/databases/uniref/uniref90_201901b_ec_filtered.fasta"
    uniref_small_output_db = "/home/chris/Projects/benchmark_HUMAnN/databases/uniref/uniref90_201901b_ec_filtered_1_small.fasta"
    uniref_large_output_db = "/home/chris/Projects/benchmark_HUMAnN/databases/uniref/uniref90_201901b_ec_filtered_2_large.fasta"


    small_db_records = (r for r in SeqIO.parse(uniref_input_db, "fasta") if r.id.split("|")[0] in small_db_whitelist)
    count = SeqIO.write(small_db_records, uniref_small_output_db, "fasta-2line")
    print(f"Wrote {count} sequences to small db")

    large_db_records = (r for r in SeqIO.parse(uniref_input_db, "fasta") if r.id.split("|")[0] not in small_db_whitelist)
    count=SeqIO.write(large_db_records, uniref_large_output_db, "fasta-2line")
    print(f"Wrote {count} sequences to large db")