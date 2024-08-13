from pathlib import Path

BENCHMARK_DIRECTORY = Path(__file__).parent / "benchmark"


def get_benchmark_directory(threads, sequences):
    return BENCHMARK_DIRECTORY / f"benchmark_{threads}threads_{sequences}seqs"


def read_benchmark_log_file(threads, sequences):
    benchmark_log_file = get_benchmark_directory(threads, sequences) / "input_data_humann_temp" / "input_data.log"
    return benchmark_log_file.read_text()
