#!/usr/bin/env python3
import multiprocessing
import queue
import shutil
import subprocess
import tarfile
from pathlib import Path
from queue import Queue
from threading import Event, Thread

import pandas as pd
from Bio import SeqIO
from manifest_processor import ManifestProcessor

files_processing_queue = Queue(maxsize=1)
download_completed_event = Event()

hmp_manifest_file = Path(__file__).parent / "hmp_manifest.tsv"
output_folder = Path(__file__).parent / "output"

default_endpoint_priority = ",".join(["HTTP", "FTP", "S3"])
diamond_database = Path(__file__).parent / "databases/uniref/uniref90_201901b_ec_filtered.dmnd"

fastq_input_relative_path = Path("input.fastq")
fasta_input_relative_path = Path("input.fasta")

cleanup = True

def extract_fastq(download_path: Path, item_folder: Path):
    fastq_file_path = item_folder / fastq_input_relative_path

    if fastq_file_path.exists():
        print("FastQ file already exists")
        print("Skipping fastq extraction")
        return fastq_file_path

    temp_fastq = fastq_file_path.with_suffix(".temp")
    with temp_fastq.open("wb") as fastq_file:
        def _recursive_extract(handle, filename):
            print("Extracting", filename)
            if filename.endswith(".fastq") or filename.endswith(".fq"):
                shutil.copyfileobj(handle, fastq_file)
            elif filename.endswith(".tar.bz2"):
                with tarfile.open(fileobj=handle, mode="r:bz2") as tar:
                    for member in tar:
                        _recursive_extract(tar.extractfile(member), member.name)
            elif filename.endswith(".tar.xz"):
                with tarfile.open(fileobj=handle, mode="r:xz") as tar:
                    for member in tar:
                        _recursive_extract(tar.extractfile(member), member.name)
            elif filename.endswith(".tar.gz"):
                with tarfile.open(fileobj=handle, mode="r:gz") as tar:
                    for member in tar:
                        _recursive_extract(tar.extractfile(member), member.name)
            elif filename.endswith(".tar"):
                with tarfile.open(fileobj=handle, mode="r") as tar:
                    for member in tar:
                        _recursive_extract(tar.extractfile(member), member.name)

        for path in download_path.iterdir():
            with path.open("rb") as handle:
                _recursive_extract(handle, path.name)

    temp_fastq.replace(fastq_file_path)
    return fastq_file_path


def extract_fasta(fastq_path: Path, item_folder: Path):
    fasta_file_path = item_folder / fasta_input_relative_path
    if fasta_file_path.exists():
        print("Fasta file already exists")
        print("Skipping fasta extraction")
        return fasta_file_path
    temp_fasta = fasta_file_path.with_suffix(".temp")
    SeqIO.convert(fastq_path, "fastq", temp_fasta, "fasta")
    temp_fasta.replace(fasta_file_path)
    return fasta_file_path


def download_files():
    mp = ManifestProcessor()
    mp.aws_s3.connection.host = "s3-accelerate.amazonaws.com"
    manifest = pd.read_csv(hmp_manifest_file, sep="\t")

    for i, single_file_manifest in manifest.iterrows():
        single_file_manifest = [{
            "id": single_file_manifest["file_id"],
            "md5": single_file_manifest["md5"],
            "urls": single_file_manifest["urls"],
        }]
        item_folder = output_folder / str(i)
        download_folder = item_folder / "download"
        fastq_file_path = item_folder / fastq_input_relative_path
        fasta_file_path = item_folder / fasta_input_relative_path
        humann_done = (item_folder / "humann").exists()
        diamond_done = (item_folder / "diamond").exists()
        fastq_exists = fastq_file_path.exists()
        fasta_exists = fasta_file_path.exists()

        if (humann_done or fastq_exists) and \
           (diamond_done or fasta_exists):
            print("Download not necessary, all required files present")
        else:
            download_folder.mkdir(parents=True, exist_ok=True)

            failed_files = mp.download_manifest(single_file_manifest, str(download_folder), default_endpoint_priority)

            if failed_files[0] != 0:
                print("Download failed!")
                if cleanup:
                    print(f"Removing {item_folder}")
                    shutil.rmtree(item_folder)
                continue

            print("Extracting fastq and fasta files...")
            fastq_file_path = extract_fastq(download_folder, item_folder)
            fasta_file_path = extract_fasta(fastq_file_path, item_folder)

        files_processing_queue.put((item_folder, fastq_file_path, fasta_file_path))
        print(f"Files for {i} enqueued")

        if cleanup:
            print(f"Removing download folder")
            shutil.rmtree(download_folder, ignore_errors=True)

    download_completed_event.set()


def run_diamond(item_folder, fasta_path):
    try:
        output_dir = item_folder / "diamond"
        if output_dir.exists():
            print("Diamond already ran!")
        else:
            temp_output_dir = output_dir.with_suffix(".temp")
            temp_output_dir.mkdir(parents=True, exist_ok=True)
            output_file = temp_output_dir / "diamond_output.tsv"

            cmd = ["diamond", "blastx",
                   "--query", str(fasta_path.absolute()),
                   "--evalue", "1.0",
                   "--threads", str(multiprocessing.cpu_count()),
                   "--top", "1",
                   "--outfmt", "6",
                   "--db", str(diamond_database.absolute()),
                   "--out", str(output_file.absolute()),
                   ]
            subprocess.check_call(cmd)
            temp_output_dir.replace(output_dir)
            print("Diamond done!")
        if cleanup:
            print("Removing input fasta")
            fasta_path.unlink()
    except:
        print(f"Diamond failed to run for folder {output_file.parent}")


def run_humann(item_folder:Path, fastq_path:Path):
    try:
        output_dir = item_folder / "humann"
        if output_dir.exists():
            print("HUMAnN already ran!")
        else:
            temp_output_dir = output_dir.with_suffix(".temp")
            temp_output_dir.mkdir(parents=True, exist_ok=True)
            cmd = ["humann",
                   "-i", str(fastq_path.absolute()),
                   "-o", str(temp_output_dir.absolute()),
                   "--resume",
                   "--threads", str(multiprocessing.cpu_count()),
#                   "--remove-temp-output",
                   ]
            subprocess.check_call(cmd)
            temp_output_dir.replace(output_dir)
            print("HUMAnN done!")

        if cleanup:
            print("Removing input fastq")
            fastq_path.unlink()
            print("Removing temp humann files")
            shutil.rmtree(output_dir / "input_humann_temp", ignore_errors=True)
    except:
        print(f"HUMAnN failed to run for folder {output_dir}")

def run_processing():
    while not download_completed_event.is_set():
        try:
            item_folder, fastq_path, fasta_path = files_processing_queue.get(timeout=5)
            run_humann(item_folder, fastq_path)
            run_diamond(item_folder, fasta_path)
        except queue.Empty:
            print("Download not ready")


if __name__ == "__main__":
    #    shutil.rmtree(output_folder, ignore_errors=True)
    output_folder.mkdir(exist_ok=True)

    process_thread = Thread(target=run_processing, daemon=False)
    process_thread.start()
    # download_thread = Thread(target=download_files)
    # download_thread.start()
    # download_thread.join()

    download_files()

    process_thread.join()

