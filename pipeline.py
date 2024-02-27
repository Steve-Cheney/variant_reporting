import subprocess
import argparse
import os
import time

def run_demultiplex(fastq, clinical_data):
    # Command to run demultiplex.py
    command = ["python3", "demultiplex.py", "-f", fastq, "-d", clinical_data]
    subprocess.run(command, check=True)

def run_alignment(reference_fasta, clinical_data):
    # Command to run align.py
    command = ["python3", "align.py", "-r", reference_fasta, "-f", "fastqs", "-d", clinical_data]
    subprocess.run(command, check=True)

def main():
    start_time = time.perf_counter()
    # Define arguments
    parser = argparse.ArgumentParser(description="Pipeline runner")
    parser.add_argument("-f", "--fastq", required=True, help="Path to FASTQ file")
    parser.add_argument("-d", "--data", required=True, help="Path to clinical data text file")
    parser.add_argument("-r", "--ref", required=True, help="Path to reference FASTA")
    args = parser.parse_args()

    # Run demultiplex.py
    run_demultiplex(args.fastq, args.data)

    # Run align.py
    run_alignment(args.ref, args.data)
    end_time = time.perf_counter()
    print(f"\n\nProcess completed in {round(end_time-start_time, 3)} seconds.")

if __name__ == "__main__":
    main()
