import argparse
import multiprocessing
import re
import shutil
import tarfile
import traceback
from collections import Counter, defaultdict
from io import TextIOWrapper
from pathlib import Path


input_folder = Path("/home/chris/chris_humann_processed_sra")

accessions = [
    "SRR15372865", "SRR15372866", "SRR15372867", "SRR15372868", "SRR15372869", "SRR15372870", "SRR15372871",
    "SRR15372873", "SRR15372874", "SRR15372875", "SRR15372876", "SRR15372877", "SRR15372878", "SRR15372879",
    "SRR15372880", "SRR15372881", "SRR15372882", "SRR15372883", "SRR15372884", "SRR15372885", "SRR15372886",
    "SRR15372887", "SRR15372888", "SRR15372889", "SRR15372890", "SRR15372891", "SRR15372892", "SRR15372893",
    "SRR15372894", "SRR15372895", "SRR15372896", "SRR15372897", "SRR15372898", "SRR15372899", "SRR15372900",
    "SRR15372901", "SRR15372902", "SRR15372903", "SRR15372904", "SRR15372905", "SRR15372906", "SRR15372907",
    "SRR15372908", "SRR15372909", "SRR15372910", "SRR15372911", "SRR15372912", "SRR15372913", "SRR15372914",
    "SRR15372915", "SRR15372916", "SRR15372917", "SRR15372918", "SRR15372919", "SRR15372920", "SRR15372921",
    "SRR15372922", "SRR15372923", "SRR15372924", "SRR15372925", "SRR15372926", "SRR15372927", "SRR15372928",
    "SRR15372929", "SRR15372930", "SRR15372931", "SRR15372932", "SRR15372933", "SRR15372934", "SRR15372935",
    "SRR15372936", "SRR15372937", "SRR15372938", "SRR15372939", "SRR15372940", "SRR15372941", "SRR15372942",
    "SRR15372943", "SRR15372944", "SRR15372945", "SRR15372946", "SRR15372947", "SRR15372948", "SRR15372949",
    "SRR15372950", "SRR15372951", "SRR15372952", "SRR15372953", "SRR15372954", "SRR15372955", "SRR15372956",
    "SRR15372957", "SRR15372958", "SRR15372959", "SRR15372960", "SRR15372961", "SRR15372962", "SRR15372963",
    "SRR15372964", "SRR15372965", "SRR15372966", "SRR15372967", "SRR15372968", "SRR15372969", "SRR15372970",
    "SRR15372971", "SRR15372972", "SRR15372973", "SRR15372974", "SRR15372975", "SRR15372976", "SRR15372977",
    "SRR15372979", "SRR15372980", "SRR15372981", "SRR15372982", "SRR15372983", "SRR15372984", "SRR15372985",
    "SRR15372986", "SRR15372987", "SRR15372988", "SRR15372989", "SRR15372990", "SRR15372991", "SRR15372992",
    "SRR15372993", "SRR15372994", "SRR15372995", "SRR15372996", "SRR15372997", "SRR15372998", "SRR15372999",
    "SRR15373000", "SRR15373001", "SRR15373002", "SRR15373003", "SRR15373004", "SRR15373005", "SRR15373006",
    "SRR15373007", "SRR15373008", "SRR15373009", "SRR15373010", "SRR15373011", "SRR15373012", "SRR15373013",
    "SRR15373014", "SRR15373015", "SRR15373016", "SRR15373017", "SRR15373018", "SRR15373019", "SRR15373020",
    "SRR15373021", "SRR15373022", "SRR15373023", "SRR15373024", "SRR15373025", "SRR15373026", "SRR15373027",
    "SRR15373028", "SRR15373029", "SRR15373030", "SRR15373031", "SRR15373032", "SRR15373033", "SRR15373034",
    "SRR15373035", "SRR15373036", "SRR15373037", "SRR15373038", "SRR15373039", "SRR15373040", "SRR15373041",
    "SRR15373042", "SRR15373043", "SRR15373044", "SRR15373045", "SRR15373046", "SRR15373047", "SRR15373048",
    "SRR15373049", "SRR15373050", "SRR15373051", "SRR15373052", "SRR15373053", "SRR15373054", "SRR15373055",
    "SRR15373056", "SRR15373057", "SRR15373058", "SRR15373059", "SRR15373060", "SRR15373061", "SRR15373062",
    "SRR15373063", "SRR15373064", "SRR15373065", "SRR15373066", "SRR15373067", "SRR15373068", "SRR15373069",
    "SRR15373070", "SRR15373071", "SRR15373072", "SRR15373073", "SRR15373074", "SRR15373075", "SRR15373076",
    "SRR15373077", "SRR15373078", "SRR15373079", "SRR15373080", "SRR15373081", "SRR15373082", "SRR15373083",
    "SRR15373084", "SRR15373085", "SRR15373086", "SRR15373087", "SRR15373088", "SRR15373089", "SRR15373090",
    "SRR15373091", "SRR15373092", "SRR15373093", "SRR15373094", "SRR15373095", "SRR15373096", "SRR15373097",
    "SRR15373098", "SRR15373099", "SRR15373100", "SRR15373101", "SRR15373102", "SRR15373103", "SRR15373104",
    "SRR15373105", "SRR15373106", "SRR15373107", "SRR15373108", "SRR15373109", "SRR15373110", "SRR15373111",
    "SRR15373112", "SRR15373113", "SRR15373114", "SRR15373115", "SRR15373116", "SRR15373117", "SRR15373118",
    "SRR15373119", "SRR15373120", "SRR15373121", "SRR15373122", "SRR15373123", "SRR15373124", "SRR15373125",
    "SRR15373126", "SRR15373127", "SRR15373128", "SRR15373129", "SRR15373130", "SRR15373131", "SRR15373132",
    "SRR15373133", "SRR15373134", "SRR15373135", "SRR15373136", "SRR15373137", "SRR15373138", "SRR15373139",
    "SRR15373140", "SRR15373141", "SRR15373142", "SRR15373143", "SRR15373144", "SRR15373145", "SRR15373146",
    "SRR15373147", "SRR15373148", "SRR15373149", "SRR15373150", "SRR15373151", "SRR15373152", "SRR15373153",
    "SRR15373154", "SRR15373155", "SRR15373156", "SRR15373157", "SRR15373158", "SRR15373159", "SRR15373160",
    "SRR15373161", "SRR15373162", "SRR15373163", "SRR15373164", "SRR15373165", "SRR15373166", "SRR15373167",
    "SRR15373168", "SRR15373169", "SRR15373170", "SRR15373171", "SRR15373172", "SRR15373173", "SRR15373174",
    "SRR15373175", "SRR15373176", "SRR15373177", "SRR15373178", "SRR15373179", "SRR15373180", "SRR15373181",
    "SRR15373182", "SRR15373183", "SRR15373184", "SRR15373185", "SRR15373186", "SRR15373187", "SRR15373188",
    "SRR15373189", "SRR15373190", "SRR15373191", "SRR15373192", "SRR15373193", "SRR15373194", "SRR15373195",
    "SRR15373196", "SRR15373197", "SRR15373198", "SRR15373199", "SRR15373200", "SRR15373201", "SRR15373202",
    "SRR23829177", "SRR23829178", "SRR23829181", "SRR23829182", "SRR23829185", "SRR23829186", "SRR23829188",
    "SRR23829191", "SRR23829202", "SRR23829205", "SRR23829208", "SRR23829210", "SRR23829211", "SRR23829212",
    "SRR23829213", "SRR23829214", "SRR23829215",
]

output_folder = Path(__file__).parent.resolve() / 'output'


def process_log(log_stream):
    log = log_stream.read()

    # b'200000 reads; of these:
    total_match = re.search(r"(\d+) reads; of these:", log)
    total_reads = int(total_match.group(1))

    # Unaligned reads after nucleotide alignment: 74.7759751156 %
    nucl_unaligned_percentage_match = re.search(r"Unaligned reads after nucleotide alignment: (\d+\.\d+) %", log)
    nucleotide_unaligned_reads = int(total_reads * float(nucl_unaligned_percentage_match.group(1)) / 100)

    # Unaligned reads after translated alignment: 98.9000000000 %
    transl_unaligned_percentage_match = re.search(r"Unaligned reads after translated alignment: (\d+\.\d+) %", log)
    translated_unaligned_reads = int(total_reads * float(transl_unaligned_percentage_match.group(1)) / 100)

    return total_reads, nucleotide_unaligned_reads, translated_unaligned_reads


def process_diamond_file(diamond_stream):
    gene_families = Counter(line.split("\t")[1].split("|")[0] for line in diamond_stream)
    return gene_families


def extract_counts(results_folder: Path):
    print(f"Parsing data for accession {results_folder.name}")
    item_output_folder = output_folder / results_folder.name
    try:
        metadata_archive_path = results_folder / "metadata.tar.xz"
        with tarfile.open(metadata_archive_path) as metadata_tar:
            log_reader = TextIOWrapper(metadata_tar.extractfile("input.log"))
            diamond_reader = TextIOWrapper(metadata_tar.extractfile("diamond_aligned.tsv"))
            total_reads, nucleotide_unaligned_reads, translated_unaligned_reads = process_log(log_reader)
            gene_family_counter = process_diamond_file(diamond_reader)

        item_output_folder.mkdir(exist_ok=True, parents=True)

        diamond_condensed_path = item_output_folder / "diamond_condensed.tsv"
        with diamond_condensed_path.open("w") as diamond_output:
            lines_iterator = (f"{gene_family}\t{gf_count}\n" for gene_family, gf_count in gene_family_counter.items())
            diamond_output.writelines(lines_iterator)

        counts_path = item_output_folder / "counts.tsv"
        with counts_path.open("w") as counts_output:
            lines = ["type\tcount\n",
                     f"total\t{total_reads}\n",
                     f"nucleotide_unaligned\t{nucleotide_unaligned_reads}\n",
                     f"translated_unaligned\t{translated_unaligned_reads}\n", ]
            counts_output.writelines(lines)

        return gene_family_counter, (total_reads, nucleotide_unaligned_reads, translated_unaligned_reads)

    except Exception as e:
        print(f"Failed for folder {results_folder}")
        traceback.print_exc()
        shutil.rmtree(item_output_folder, ignore_errors=True)
        return dict(), (0,0,0)


def main():
    total_gene_families_count = defaultdict(int)
    total_reads_count = 0
    total_nucleotide_unaligned_count = 0
    total_translated_unaligned_count = 0

    pool = multiprocessing.Pool()

    for gene_family_counter, (total_reads, nucleotide_unaligned_reads, translated_unaligned_reads) \
        in pool.imap_unordered(extract_counts, (input_folder / accession for accession in accessions)):

        for gene_family, count in gene_family_counter.items():
            total_gene_families_count[gene_family] += count

        total_reads_count += total_reads
        total_nucleotide_unaligned_count += nucleotide_unaligned_reads
        total_translated_unaligned_count += translated_unaligned_reads

    total_gene_families_count = list(total_gene_families_count.items())
    total_gene_families_count.sort(reverse=True, key=lambda x: x[1])

    with (output_folder / "total_aa_counts.tsv").open("w") as f:
        f.write("gene_family\tcount\n")
        for amino_acid, count in total_gene_families_count:
            f.write(f"{amino_acid}\t{count}\n")

    with (output_folder / "unaligned_counts.tsv").open("w") as f:
        lines = ["type\tcount\n",
                 f"total\t{total_reads_count}\n",
                 f"nucleotide_unaligned\t{total_nucleotide_unaligned_count}\n",
                 f"translated_unaligned\t{total_translated_unaligned_count}\n", ]
        f.writelines(lines)


if __name__ == '__main__':
    main()

    # with open("/home/chris/Projects/benchmark_HUMAnN/benchmark_uniref90_full/benchmark_20threads_200000seqs/input_data_humann_temp/input_data.log") as input_handle:
    #    print(process_log(input_handle))

    # with open("/home/chris/Projects/benchmark_HUMAnN/benchmark_uniref90_full/benchmark_20threads_200000seqs/input_data_humann_temp/input_data_diamond_aligned.tsv") as input_handle:
    #    print(process_diamond_file(input_handle))
