from Bio import SeqIO
import os

read1_file = "/home/roshanpe/scGESTALT/GESTALT_FINAL_SCRIPT/Brain_1_1597_1083_polyA_filtered_reads/polyA_filtered_read1.fastq"
read2_file = "/home/roshanpe/scGESTALT/GESTALT_FINAL_SCRIPT/Brain_1_1597_1083_polyA_filtered_reads/polyA_filtered_read2.fastq"
output_folder = "polyA_1597_1083_after_short_read2_filtered_reads"

# Create the output folder if it doesn't exist
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Read read1 sequences into a dictionary
read1_records_dict = {record.id: record for record in SeqIO.parse(read1_file, "fastq")}

# Read read2 sequences into a list and filter out sequences with length < 10
read2_records = [record for record in SeqIO.parse(read2_file, "fastq") if len(record.seq) >= 10]

# Filter read1 sequences to match read2 and remove those with length < 10
read1_records_filtered = [read1_records_dict[record.id] for record in read2_records]

# Write all filtered read1 sequences to a single fastq file
output_file_read1 = os.path.join(output_folder, "polyA_after_short_reads_filter_read1.fastq")
with open(output_file_read1, "w") as output_handle_read1:
    SeqIO.write(read1_records_filtered, output_handle_read1, "fastq")

# Write all filtered read2 sequences to a single fastq file
output_file_read2 = os.path.join(output_folder, "polyA_after_short_reads_filter_read2.fastq")
with open(output_file_read2, "w") as output_handle_read2:
    SeqIO.write(read2_records, output_handle_read2, "fastq")

