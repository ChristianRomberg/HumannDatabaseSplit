#!/usr/bin/env python3
import itertools
import random
import shutil
import subprocess
from pathlib import Path

from Bio import SeqIO

from util import get_benchmark_directory, read_benchmark_log_file

NUM_THREADS_TO_TEST = [20]
NUM_SEQUENCES_TO_TEST = [1000, 20000, 50000, 100000, 200000, 500000, 1000000, 2000000, 5000000]
# INPUT_FILE_PATH = Path(__file__).parent / "input_data" / "humann_benchmark.fasta"
INPUT_FILE_PATH = Path("input_data/CL100028205_L02_40.rmhost_kneaddata.fastq")
ALL_SEQUENCES_INDEX = SeqIO.index(str(INPUT_FILE_PATH), "fastq")
ALL_SEQUENCES_KEYS = list(ALL_SEQUENCES_INDEX.keys())

RANDOM_SEED = 42


def run_benchmark(num_threads_to_test, num_sequences_to_test):
    try:
        benchmark_log = read_benchmark_log_file(num_threads_to_test, num_sequences_to_test)
        assert "Output files created" in benchmark_log
        print(f"Benchmark for {num_threads_to_test} threads and {num_sequences_to_test} sequences already run!")
        return
    except (AssertionError, FileNotFoundError):
        # Benchmark not yet run
        pass

    print(f"Preparing benchmark for {num_threads_to_test} threads and {num_sequences_to_test} sequences")

    output_dir = get_benchmark_directory(num_threads_to_test, num_sequences_to_test)
    shutil.rmtree(output_dir, ignore_errors=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    random.seed(RANDOM_SEED)
    selected_sequences = random.sample(ALL_SEQUENCES_KEYS, k=num_sequences_to_test)
    input_file = output_dir / "input_data.fastq"
    SeqIO.write((ALL_SEQUENCES_INDEX[key] for key in selected_sequences),
                input_file, "fastq")

    print(f"Benchmark starting")

    cmd = ["humann",
           "-i", str(input_file.absolute()),
           "-o", str(output_dir.absolute()),
           "--threads", str(num_threads_to_test)]

    subprocess.check_call(cmd)


def compute_average_seq_length(average_seq_length_number=10000):
    selected_sequences = random.sample(ALL_SEQUENCES_KEYS, k=average_seq_length_number)
    len_sum = sum(len(ALL_SEQUENCES_INDEX[seq_key].seq) for seq_key in selected_sequences)
    return len_sum / average_seq_length_number


def main():
    print(f"Average sequence length: {compute_average_seq_length()}")
    for num_threads_to_test, num_sequences_to_test in itertools.product(NUM_THREADS_TO_TEST, NUM_SEQUENCES_TO_TEST):
        run_benchmark(num_threads_to_test, num_sequences_to_test)


if __name__ == '__main__':
    main()
