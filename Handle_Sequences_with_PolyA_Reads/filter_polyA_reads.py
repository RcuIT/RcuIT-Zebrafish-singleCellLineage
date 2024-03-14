from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os

# Function to find poly-A sequences and trim them
def trim_poly_a(sequence):
    index = sequence.find("AAAAA")  # Find the index of the first "AAAAA" sequence
    if index != -1:
        return sequence[:index]  # Trim the sequence till the first "A" of "AAAAA"
    else:
        return sequence  # Return the original sequence if "AAAAA" is not found

# Paths to input and output files
read1_file = "/home/roshanpe/scGESTALT/GESTALT_FINAL_SCRIPT/Brain.1.1597_1083_unommon_CB/Brain.1.1597_1083_unommon_read1.trimmed.fastq"
read2_file = "/home/roshanpe/scGESTALT/GESTALT_FINAL_SCRIPT/Brain.1.1597_1083_unommon_CB/Brain.1.1597_1083_unommon_read2.filtered.fastq"
output_folder = "Brain_1_1597_1083_polyA_filtered_reads"

# Create the output folder if it doesn't exist
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Lists to store trimmed read1 and read2 sequences
trimmed_read1_records = []
trimmed_read2_records = []

# Open read2 file for parsing and writing trimmed reads
with open(read2_file, "r") as read2_handle, open(read1_file, "r") as read1_handle:
    read1_dict = SeqIO.to_dict(SeqIO.parse(read1_handle, "fastq"))
    for record in SeqIO.parse(read2_handle, "fastq"):
        trimmed_sequence = trim_poly_a(str(record.seq))
        if trimmed_sequence != str(record.seq):  # Check if trimming occurred
            # Create a SeqRecord object for trimmed read2 sequence
            read2_record = SeqRecord(Seq(trimmed_sequence), id=record.id, description="")
            # Assign a default quality score of 40 to the trimmed sequence
            read2_record.letter_annotations["phred_quality"] = [40] * len(trimmed_sequence)
            trimmed_read2_records.append(read2_record)
            # Get the header from read2
            read1_header = record.id.split()[0]
            trimmed_read1_records.append(read1_dict[read1_header])

# Write all trimmed read1 sequences to one file
output_file_read1 = os.path.join(output_folder, "polyA_filtered_read1.fastq")
with open(output_file_read1, "w") as output_handle:
    SeqIO.write(trimmed_read1_records, output_handle, "fastq")

# Write all trimmed read2 sequences to another file
output_file_read2 = os.path.join(output_folder, "polyA_filtered_read2.fastq")
with open(output_file_read2, "w") as output_handle:
    SeqIO.write(trimmed_read2_records, output_handle, "fastq")