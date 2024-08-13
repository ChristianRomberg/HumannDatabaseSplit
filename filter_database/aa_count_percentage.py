#!/usr/bin/env python3
import bisect

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import StrMethodFormatter

if __name__ == '__main__':
    total_gene_families = 6191048
    aa_counts = pd.read_csv("total_aa_counts.tsv", sep="\t", index_col=0, header=None, names=["protein", "count"])

    total = aa_counts["count"].sum()

    cumsum_fraction = np.array(aa_counts["count"].cumsum() / total * 100)
    db_size_fraction = (np.arange(len(cumsum_fraction)) + 1) / total_gene_families

    time_first_stage = db_size_fraction * 100
    time_second_stage = (1 - (cumsum_fraction / 100)) * (1 - db_size_fraction) * 100
    # time_fraction = 100 * ((db_size_fraction) + (1 - (cumsum_fraction / 100)) * (1 - db_size_fraction))
    time_both_stages = time_first_stage + time_second_stage
    quantiles_y = [i * 10 for i in range(1, 11)]
    quantiles_x = [bisect.bisect_left(cumsum_fraction, quantile) for quantile in quantiles_y]

    plt.figure(figsize=(8.4, 4.6), dpi=72 * 2)
    # plt.semilogx(cumsum_fraction)
    plt.plot(cumsum_fraction, label="% reads covered")
    plt.plot(time_first_stage, label="% time for first stage")
    plt.plot(time_second_stage, label="% time for second stage")
    plt.plot(time_both_stages, label="% time for both stages")
    plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))

    for quantile in quantiles_y:
        plt.axhline(quantile, linestyle=":", color="black", alpha=0.3)
    for quantile in quantiles_x:
        plt.axvline(quantile, linestyle="--", color="black", alpha=0.3)
    # plt.hlines(quantiles_y, -len(cumsum_fraction), 2*len(cumsum_fraction), color="black", alpha=0.3)
    # plt.vlines(quantiles_x, -1, 1, color="black", alpha=0.3)

    plt.xlabel("Number of gene families")
    #plt.ylabel("Percentage of mapped reads")
    #plt.title("Gene family vs cumulative abundance")
    plt.legend()
    plt.tight_layout()
    plt.show()

    print("Quantiles:")
    for quantile, gene_families in zip(quantiles_y, quantiles_x):
        print(f"{quantile}%: {gene_families}")

    print("Done!")
