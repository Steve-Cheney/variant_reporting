import argparse
import subprocess
import os
import pysam

def _pileup(sorted_bam, ref_fa) -> list:
    """
    Given a sorted bam file and reference sequence, get the SNP details.

    Params
    ------
    sorted_bam
        The {name}.sorted.bam file to pilup
    ref_fa
        The reference fasta file

    Returns
    -------
    A list of SNPs where each item in the list is a dictionary in the following format:
    [(read_count, position, ref_bp, (pileup_SNP, frequency)), ...]
    """
    ref_seq = ""
    with open(ref_fa, 'r') as f:
        # Skip the first line (header)
        next(f)
        for line in f:
            ref_seq += line.strip()
        
    #test file, replaced with the sorted.bam you are using. Make sure it is indexed! (Use samtools index yourbam.sorted.bam)
    samfile = pysam.AlignmentFile(sorted_bam, "rb")

    #Since our reference only has a single sequence, we're going to pile up ALL of the reads. Usually you would do it in a specific region (such as chromosome 1, position 1023 to 1050 for example)
    SNPs_list = []
    
    for pileupcolumn in samfile.pileup():
        #print ("coverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
        #use a dictionary to count up the bases at each position
        ntdict = {'A':0, 'C':0, 'T':0, 'G':0}
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                # You can uncomment the below line to see what is happening in the pileup. 
                #print ('\tbase in read %s = %s' % (pileupread.alignment.query_name, pileupread.alignment.query_sequence[pileupread.query_position]))
                base = pileupread.alignment.query_sequence[pileupread.query_position]
                ntdict[base] += 1
        #print (ntdict)
        freq_dict = {key: value/pileupcolumn.n for key, value in ntdict.items()}
        #print(freq_dict)
        SNP_ntdict = {key: value for key, value in freq_dict.items() if value != 1 and value != 0}
        if SNP_ntdict:
            #print("SNP found:")
            ref_bp = ref_seq[pileupcolumn.pos]
            SNP_tuple = (None, None)
            for key, value in SNP_ntdict.items():
                if key != ref_bp:
                    SNP_tuple = (key, value)
            SNP_details = pileupcolumn.n, pileupcolumn.pos, ref_bp, SNP_tuple
            SNPs_list.append(SNP_details)
    samfile.close()
    return SNPs_list

def run_alignment(folder_path, ref_path):
    """
    Align the sam files and create the appropriate bam files for variant analysis.

    Params
    ------
    folder_path
        The path to the fastqs folder
    ref_path
        The path to the reference fasta file

    Returns
    -------
    A list of of all details of the alignments and SNPs
    """
    # Check if the folder exists
    if not os.path.exists(folder_path):
        print("Folder does not exist.")
        return
    
    # Check if the reference file exists
    if not os.path.exists(ref_path):
        print("Reference file does not exist.")
        return
    
    # Create a subfolder named "sam" within the given folder_path
    sam_folder = os.path.join(folder_path, 'sam')
    os.makedirs(sam_folder, exist_ok=True)
    
    # Create "bam" folder in the same directory as fastqs
    parent_dir = os.path.dirname(folder_path)
    bam_folder = os.path.join(parent_dir, 'bam')
    os.makedirs(bam_folder, exist_ok=True)
    
    output = []

    # Iterate over files in the folder
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        # Check if the path is a file
        if os.path.isfile(file_path):
            file_name_parts = filename.split("_trimmed.fastq")
            if len(file_name_parts) == 2:
                name = file_name_parts[0]
                # Run BWA alignment command
                sam_file_path = os.path.join(sam_folder, f'{name}.sam')
                with open(sam_file_path, 'w') as output_file:
                    subprocess.run(['bwa', 'mem', ref_path, file_path], stdout=output_file)
                
                # Convert SAM to BAM
                bam_file_path = os.path.join(bam_folder, f'{name}.bam')
                with open(bam_file_path, 'wb') as output_file:
                    subprocess.run(['samtools', 'view', '-bS', '-o', '-', sam_file_path], stdout=output_file)

                # Sort BAM file
                sorted_bam_file_path = os.path.join(bam_folder, f'{name}.sorted.bam')
                subprocess.run(['samtools', 'sort', '-m', '100M', '-o', sorted_bam_file_path, bam_file_path])

                # Index sorted BAM file
                subprocess.run(['samtools', 'index', sorted_bam_file_path])
                
                snps = _pileup(sorted_bam_file_path, ref_path)
                output.append({name:snps})
    return output
                
def cd_to_dict(file_path):
    """
    Given a clinical data text file, convert it to a dict.

    Params
    ------
    file_path
        The path to the clinical data file
   
    Returns
    -------
    A dictionary in the format {name: color}
    """    
    data_dict = {}
    with open(file_path, 'r') as file:
        # Skip the header line
        next(file)
        for line in file:
            name, color, barcode = line.strip().split('\t')
            data_dict[name] = color
    return data_dict

def write_report(sample_details, output_file):
    with open(output_file, 'w') as f:
        for sample in sample_details:
            name, details = list(sample.items())[0]
            mutation_info = details[0]
            color = details[1]
            reads = mutation_info[0]
            position = mutation_info[1]+1
            mutation = mutation_info[3][0]
            mutation_percentage = mutation_info[3][1] * 100

            f.write(f"Sample {name} had a {color} mold, {reads} reads, ")
            f.write(f"and {mutation_percentage:.1f}% of the reads at position {position} ")
            f.write(f"had the mutation {mutation}.\n")

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--ref", required=True, help="Place path to reference FASTA inside here")
    parser.add_argument("-f", "--folder", required=True, help="Place path to FASTQ folder inside here")
    parser.add_argument("-d", "--data", required=True, help="Place clinical data text file inside here")
    
    args = parser.parse_args()

    ref_fa = args.ref
    fastq_folder = args.folder

    clinical_data_file = args.data
    cd_dict = cd_to_dict(clinical_data_file)

    subprocess.run(['bwa', 'index', ref_fa])
    details = run_alignment(fastq_folder, ref_fa)
    for detail in details:
        for key, value in detail.items():
            detail[key].append(cd_dict[key])

    write_report(details, "report.txt")
