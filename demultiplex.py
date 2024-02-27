import argparse
import os
import shutil

class ParseFastQ(object):
    """Returns a read-by-read fastQ parser analogous to file.readline()"""
    def __init__(self,filePath,headerSymbols=['@','+']):
        """Returns a read-by-read fastQ parser analogous to file.readline().
        Exmpl: parser.next()
        -OR-
        Its an iterator so you can do:
        for rec in parser:
            ... do something with rec ...
 
        rec is tuple: (seqHeader,seqStr,qualHeader,qualStr)
        """
        if filePath.endswith('.gz'):
            self._file = gzip.open(filePath)
        else:
            self._file = open(filePath, 'r')
        self._currentLineNumber = 0
        self._hdSyms = headerSymbols
         
    def __iter__(self):
        return self
     
    def __next__(self):
        """Reads in next element, parses, and does minimal verification.
        Returns: tuple: (seqHeader,seqStr,qualHeader,qualStr)"""
        # ++++ Get Next Four Lines ++++
        elemList = []
        for i in range(4):
            line = self._file.readline()
            self._currentLineNumber += 1 ## increment file position
            if line:
                elemList.append(line.strip('\n'))
            else: 
                elemList.append(None)
         
        # ++++ Check Lines For Expected Form ++++
        trues = [bool(x) for x in elemList].count(True)
        nones = elemList.count(None)
        # -- Check for acceptable end of file --
        if nones == 4:
            raise StopIteration
        # -- Make sure we got 4 full lines of data --
        assert trues == 4,\
               "** ERROR: It looks like I encountered a premature EOF or empty line.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber)
        # -- Make sure we are in the correct "register" --
        assert elemList[0].startswith(self._hdSyms[0]),\
               "** ERROR: The 1st line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[0],self._currentLineNumber) 
        assert elemList[2].startswith(self._hdSyms[1]),\
               "** ERROR: The 3rd line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[1],self._currentLineNumber) 
        # -- Make sure the seq line and qual line have equal lengths --
        assert len(elemList[1]) == len(elemList[3]), "** ERROR: The length of Sequence data and Quality data of the last record aren't equal.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber) 
         
        # ++++ Return fatsQ data as tuple ++++
        return tuple(elemList)

def cd_to_dict(file_path):
    """
    Given a clinical data text file, convert it to a dict.

    Params
    ------
    file_path
        The path to the clinical data file
   
    Returns
    -------
    A dictionary in the format {barcode: 'Name': name, 'Color': color}
    """   

    data_dict = {}
    with open(file_path, 'r') as file:
        # Skip the header line
        next(file)
        for line in file:
            name, color, barcode = line.strip().split('\t')
            data_dict[barcode] = {'Name': name, 'Color': color}
    return data_dict

def trim_reads(seq, quality_score):
    """
    Trim the reads of a sequence and qs if the end of the read has at least two consecutive quality scores of D or F.

    Params
    ------
    seq
        The sequence to trim
    quality_score
        The quality score to reference and also trim
   
    Returns
    -------
    The trimmed tuple of sequence and quality score
    """     

    qs_length = len(quality_score) 
    for i, bp in enumerate(quality_score):
        if bp == 'D' or bp == 'F':
            if i == qs_length - 1:
                return seq, quality_score
            next_char = quality_score[i+1] 
            if next_char == 'D' or next_char == 'F':
                return seq[0:i], quality_score[0:i]

def append_fastq_to_file(cd_dict, barcode, string_to_append):
    """
    Append a given fastq sequence to the appropriate file.

    Params
    ------
    cd_dict
        The clinical data dictionary from cd_to_dict()
    barcode
        The barcode given from the clinical data text
    string_to_append
        The fastq sequence to append

    Returns
    -------
    
    """   

    subfolder = "fastqs"
    if barcode in cd_dict:
        name = cd_dict[barcode]["Name"]
        filename = os.path.join(subfolder, f"{name}_trimmed.fastq")
        os.makedirs(subfolder, exist_ok=True)
        with open(filename, "a+") as file:
            file.write(string_to_append)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fastq", required=True, help="Place fastq inside here")
    parser.add_argument("-d", "--data", required=True, help="Place clinical data text file inside here")
    args = parser.parse_args()

    # Parse the args 
    fastqfile = ParseFastQ(args.fastq)
    clinical_data_file = args.data

    cd_dict = cd_to_dict(clinical_data_file)
    #print(cd_dict)
    

    # Header: fastq_obj[0]
    # Sequence: fastq_obj[1]
    # Separator: fastq_obj[2]
    # Quality score: fastq_obj[3]
    try:
        # Attempt to remove any existing fastqs directory and all its contents
        shutil.rmtree('fastqs')
        print(f"Existing directory '{'fastqs'}' and its contents have been removed.")
    except Exception:
        pass

    for fastq_obj in fastqfile:
        header = fastq_obj[0]
        seq = fastq_obj[1]
        separator = fastq_obj[2]
        q_score = fastq_obj[3]

        barcode = seq[0:5]
        headless_seq = seq[5:]
        headless_q_score = q_score[5:]
        
        trimmed_seq, trimmed_q_score = trim_reads(headless_seq, headless_q_score)

        output = f"{header}\n{trimmed_seq}\n{separator}\n{trimmed_q_score}\n"
        
        append_fastq_to_file(cd_dict, barcode, output)