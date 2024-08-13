import os
import pathlib
from collections import defaultdict


def count_amino_acids(files):
    total_counts = defaultdict(int)

    for file in files:
        if file.is_file() and file.suffix == '.tsv':
            print(f"Reading {file}")
            with file.open('r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    fields = line.strip().split('\t')
                    amino_acid = fields[0]
                    count = int(fields[1])
                    total_counts[amino_acid] += count
    return dict(total_counts)

if __name__ == '__main__':
    files = pathlib.Path("output").glob("*/diamond/diamond_condensed.tsv")
    total_amino_acids_count = count_amino_acids(files)
    total_amino_acids_count = list(total_amino_acids_count.items())
    total_amino_acids_count.sort(reverse=True, key=lambda x: x[1])
    with open("total_aa_counts.tsv", "w") as f:
        for amino_acid, count in total_amino_acids_count:
            f.write(f"{amino_acid}\t{count}\n")
    print("Done counting!")
    #print(total_amino_acids_count)
