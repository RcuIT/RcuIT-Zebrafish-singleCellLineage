from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os

def main():
    #text_file = "/home/roshanpe/scGESTALT/GESTALT_FINAL_SCRIPT/Brain_2_Filtered_Fq.with.Uniuq.CB.txt"
    #text_file = "/home/roshanpe/scGESTALT/GESTALT_FINAL_SCRIPT/Brain.1.1083.CB.txt"
    text_file = "/home/roshanpe/scGESTALT/GESTALT_FINAL_SCRIPT/Brain.1.1597_1083_unommon.filtered.CB.Tr.Lin.txt"
    read1_file = "/home/roshanpe/my_test_run/data/raw/Brain_1_read1.trimmed.fastq"
    read2_file = "/home/roshanpe/my_test_run/data/raw/Brain_1_read2.filtered.fastq"
    #read1_file = "/home/roshanpe/my_test_run/data/raw/Brain_1_read1.trimmed.fastq"
    #read2_file = "/home/roshanpe/my_test_run/data/raw/Brain_1_read2.filtered.fastq"

    common_output_dir = "Brain.1.1597_1083_unommon_CB"
    os.makedirs(common_output_dir, exist_ok=True)
    read1_trimmed_out_path = os.path.join(common_output_dir, "Brain.1.1597_1083_unommon_read1.trimmed.fastq")
    read2_filtered_out_path = os.path.join(common_output_dir, "Brain.1.1597_1083_unommon_read2.filtered.fastq")
    

    # Read cell barcodes from the text file
    barcode_dict = {}
    with open(text_file, "r") as txt:
        next(txt)  # Skip the header line
        for line in txt:
            parts = line.strip().split("\t")
            barcode_dict[parts[0]] = parts[1]

    # Get the set of common headers between read1 and read2
    common_headers = set(header for header in barcode_dict.keys())

    # Process read2 file, filter common sequences, and write filtered read2
    common_read2_records = set()
    with open(read2_filtered_out_path, "w") as read2_filtered_out:
        for read2_record in SeqIO.parse(read2_file, "fastq"):
            if read2_record.id in common_headers:
                SeqIO.write(read2_record, read2_filtered_out, "fastq")
                common_read2_records.add(read2_record.id)

    # Process read1 file, filter common sequences, and write replaced read1
    with open(read1_trimmed_out_path, "w") as read1_trimmed_out:
        for read1_record in SeqIO.parse(read1_file, "fastq"):
            if read1_record.id in common_read2_records:
                first_28_bases = str(read1_record.seq[:28])
                replaced_read1_seq = barcode_dict[read1_record.id] + first_28_bases[16:]
                replaced_read1_record = SeqRecord(Seq(replaced_read1_seq), id=read1_record.id, description="")
                replaced_read1_record.letter_annotations["phred_quality"] = read1_record.letter_annotations["phred_quality"][:len(replaced_read1_record.seq)]
                SeqIO.write(replaced_read1_record, read1_trimmed_out, "fastq")

if __name__ == "__main__":
    main()

