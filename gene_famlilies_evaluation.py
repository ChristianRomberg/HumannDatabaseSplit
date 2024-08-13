#!/usr/bin/env python3

if __name__ == '__main__':
    gene_families = set()
    with open("uniref90_annotated_v201901b_full/diamond_output.tsv") as diamond_output:
        for line in diamond_output:
            gene_family = line.split("\t")[1]
            gene_families.add(gene_family)

    print(len(gene_families))