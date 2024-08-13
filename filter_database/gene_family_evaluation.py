import gzip
import itertools
import random
from pathlib import Path

import matplotlib.pyplot as plt

from filter_database.hmp_gene_families_and_ec_numbers import output_folder


def load_ec_map():
    result = dict()
    i = 0
    with gzip.open(ec_map_database, "rt") as handle:
        for line in handle:
            i += 1
            ec, *uniref_numbers = line.strip().split("\t")
            for uniref_number in uniref_numbers:
                result[uniref_number] = ec
    print(f"Loaded {i} ec numbers for {len(result)} gene families")
    return result


ec_map_database = Path(__file__).parent.parent / "databases/utility_mapping/map_level4ec_uniref90.txt.gz"
ec_map = load_ec_map()


def read_gene_families(gene_family_path: Path):
    gene_families = set()

    with gene_family_path.open() as input_file:
        input_file.readline()

        for line in input_file:
            gene_family, abundance_rpk = line.strip().split('\t')
            gene_family = gene_family.split('|')[0]
            gene_families.add(gene_family)
    return gene_families


def evaluate_humann():
    gene_families_sets = []
    EC_sets = []
    gene_families_not_mapped = set()
    for folder in output_folder.iterdir():
        humann_genefamily_file = folder / "humann/input_genefamilies.tsv"
        if not humann_genefamily_file.exists():
            continue

        gene_families = read_gene_families(humann_genefamily_file)
        gene_families.remove("UNMAPPED")
        ec_families = {ec_map.get(gene_family) for gene_family in gene_families}
        gene_families_not_mapped.update({gene_family for gene_family in gene_families if gene_family not in ec_map})
        ec_families.remove(None)

        print(len(gene_families))
        print(len(ec_families))
        gene_families_sets.append(gene_families)
        EC_sets.append(ec_families)

    set_size_plot(gene_families_sets)
    set_size_plot(EC_sets)


def evaluate_diamond():
    gene_families_sets = []
    EC_sets = []
    gene_families_not_mapped = set()

    for folder in output_folder.iterdir():
        file = folder / "diamond/diamond_condensed.tsv"
        if not file.exists():
            continue

        with file.open() as input_file:
            #print(input_file.readline())
            gene_families = set(line.split("\t")[0] for line in input_file.readlines())
            ec_families = {ec_map.get(gene_family) for gene_family in gene_families}
            gene_families_not_mapped.update({gene_family for gene_family in gene_families if gene_family not in ec_map})

            print(len(gene_families))
            gene_families_sets.append(gene_families)
            EC_sets.append(ec_families)

    set_size_plot(gene_families_sets)
    set_size_plot(EC_sets)


def main():
    #evaluate_humann()
    evaluate_diamond()


def set_size_plot(sets, num_lines=100, jitter_factor=0.01):
    x = range(1, len(sets) + 1)
    sets = list(sets)
    jitter = (1 + random.gauss() * jitter_factor for _ in itertools.count())
    for i in range(num_lines):
        random.shuffle(sets)
        y = []
        total_set = set(sets[0])
        for sample in sets:
            #total_set.update(sample)
            total_set.intersection_update(sample)
            y.append(len(total_set) * next(jitter))
        plt.plot(x, y, alpha=0.2, color='blue')
    plt.xlabel("Number of samples")
    plt.ylabel("Size of total set")
    plt.show()


if __name__ == '__main__':
    main()
