**Pipeline README**

**Introduction:**

This pipeline script (`pipeline.py`) is designed to automate the process of running two Python scripts, `demultiplex.py` and `align.py`, which analyze the given data and identify the mutations that resulted in the different color molds.


**Requirements:**

To run this pipeline, you need:

1. Python 3 installed on your system.
2. Input files: FASTQ file containing sequencing data and a clinical data text file.
3. A reference FASTA file for alignment.

**Usage:**

1. **Setup:** Place `pipeline.py`, `demultiplex.py`, and `align.py` scripts in the same directory. Ensure that the required input files (FASTQ file, clinical data text file, and reference FASTA file) are also available in this directory.

2. **Running the Pipeline:** Execute the pipeline by running the `pipeline.py` script with appropriate command-line arguments:

   ```
   python3 pipeline.py -f <path_to_FASTQ_file> -d <path_to_clinical_data_file> -r <path_to_reference_FASTA_file>
   ```

   Example Usage:

   ```
   python3 pipeline.py -f hawkins_pooled_sequences.fastq -d harrington_clinical_data.txt -r dgorgon_reference.fa
   ```

3. **Description of Arguments:**
   
   - `-f` or `--fastq`: Path to the FASTQ file containing sequencing data.
   - `-d` or `--data`: Path to the clinical data text file.
   - `-r` or `--ref`: Path to the reference FASTA file for alignment.

4. **Output:** The pipeline will execute `demultiplex.py` and `align.py` scripts with the provided arguments. Output files generated by these scripts will be stored in the current directory:

   - `fastqs`: A folder containing the individual trimmed fastq files for each name within the `clinical_data_file`, created within the same directory of the script.
        - `sam`: A subfolder containing the sam files from each fastq file.
   - `bam`: A folder containing the standard, sorted, and indexed bam files fron the `sam` files generated, created within the same directory of the script.
   - `indexed reference fa` files: `bwa` indexed referene file for script execution, created within the same directory of the script.
   - `report.txt`: The final report of the variant analysis, created within the same directory of the script.

5. **External Dependencies**

    - `pysam`
    - `bwa`
    - `samtools`

**Notes:**

- Ensure that the necessary dependencies for `demultiplex.py` and `align.py` scripts are installed in your Python environment.
- The pipeline assumes that `demultiplex.py` and `align.py` scripts are correctly implemented and functional.
- Always provide correct paths to input files and reference files to ensure proper execution of the pipeline.