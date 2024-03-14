from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.Align import substitution_matrices
import csv
from collections import defaultdict

def align_sequences(reference_seq, edited_seq):
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")  # using BLOSUM62 matrix
    aligner.open_gap_score = -10  # Gap open penalty
    aligner.extend_gap_score = -0.5  # Gap extension penalty
    aligner.end_open_gap_score = -10  # End gap open penalty
    aligner.end_extend_gap_score = -0.5  # End gap extension penalty
    alignments = aligner.align(reference_seq, edited_seq)
    best_alignment = max(alignments, key=lambda x: x.score)
    ref_aligned = best_alignment[0]
    edited_aligned = best_alignment[1]

    return ref_aligned, edited_aligned

def find_differences(reference_aligned, edited_aligned):
    differences = []
    diff_start = None
    diff_type = None
    diff_count = 0
    for i, (ref_base, edited_base) in enumerate(zip(reference_aligned, edited_aligned)):
        if ref_base != edited_base:
            if edited_base == "-":
                diff_type = "D"
            else:
                diff_type = "I"
            if diff_start is None:
                diff_start = i + 1
            diff_count += 1
        elif diff_start is not None:
            differences.append((diff_type, diff_count, diff_start))
            diff_start = None
            diff_count = 0

    return differences

def format_differences(differences):
    formatted_deletions = []
    formatted_mismatches = []
    for diff_type, diff_count, diff_start in differences:
        if diff_type == "D":
            formatted_deletions.append(f"{diff_count}D+{diff_start}")
        elif diff_type == "I":
            formatted_mismatches.append(f"{diff_count}I+{diff_start}")

    return formatted_deletions, formatted_mismatches

def main():
    reference_seq_file = "reference.fasta"
    read1_seq_file = "/home/roshanpe/scGESTALT/GESTALT_FINAL_SCRIPT/linBrain_2_after_CB_correction/read1.trimmed.cb.corrected.fastq"
    read2_seq_file = "/home/roshanpe/scGESTALT/GESTALT_FINAL_SCRIPT/linBrain_2_after_CB_correction/read2.filtered.cb.corrected.fastq"

    reference_seq = SeqIO.read(reference_seq_file, "fasta").seq
    reference_seq = Seq(reference_seq.upper())

    target_sites = [
        ("GCCCTGAAGCTGAAGGACGGCGG", 49, 66, 71),
        ("GGAGTTCAAGTCCATCTACATGG", 76, 93, 98),
        ("GCTGCCCGGCTACTACTACGTGG", 103, 120, 125),
        ("GGACTGGAGGACTTCTGGGGAGG", 130, 147, 152)
    ]
    read1_records = {}
    for read1_record in SeqIO.parse(read1_seq_file, "fastq"):
        read1_records[read1_record.description.split()[0]] = read1_record

    with open("target_table.csv", "w", newline="") as target_csv_file, \
            open("collapse.csv", "w", newline="") as collapse_csv_file, \
            open("mismatch_labels.csv", "w", newline="") as mismatch_csv_file:  # Open mismatch labels CSV file
        target_csv_writer = csv.writer(target_csv_file)
        collapse_csv_writer = csv.writer(collapse_csv_file)

        target_header = ["Fastq Header", "Cell Barcode", "UMI", "Reference Sequence", "Edited Sequence", "Target 1", "Target 2", "Target 3", "Target 4"]
        collapse_header = ["Fastq Header", "Cell Barcode", "UMI", "Collapse Count", "Reference Sequence", "Edited Sequence", "Target 1", "Target 2", "Target 3", "Target 4"]

        target_csv_writer.writerow(target_header)
        collapse_csv_writer.writerow(collapse_header)

        collapsed_data = defaultdict(list)

        for read2_record in SeqIO.parse(read2_seq_file, "fastq"):
            read2_header = read2_record.description.split()[0]
            if read2_header in read1_records:
                read1_record = read1_records[read2_header]
                edited_seq = read2_record.seq
                edited_seq = Seq(edited_seq.upper())

                ref_aligned, edited_aligned = align_sequences(reference_seq, edited_seq)
                differences = find_differences(ref_aligned, edited_aligned)
                formatted_deletions, formatted_mismatches = format_differences(differences)

                target_labels = ["NONE"] * len(target_sites)

                for i, (target, start, cut_pos, end_pos) in enumerate(target_sites):
                    has_deletion = False
                    is_edited = False  # Add this variable to track if the target site is edited
                    processed_positions = set()

                    for formatted_del in formatted_deletions:
                        del_count, del_start = map(int, formatted_del.split("D+"))
                        if start <= del_start <= cut_pos:
                            if del_count >= cut_pos - start:
                                if del_count >= 81:  # checking edition among ALL adjacent targets
                                    for j in range(i, min(i + 4, len(target_labels))):
                                        target_labels[j] = formatted_del
                                        has_deletion = True
                                        processed_positions.add(del_start)
                                    break
                                elif del_count >= 54:  # checking edition among THREE adjacent targets
                                    for j in range(i, min(i + 3, len(target_labels))):
                                        target_labels[j] = formatted_del
                                        has_deletion = True
                                        processed_positions.add(del_start)
                                    break
                                elif del_count >= 27:  # checking edition among TWO adjacent targets
                                    for j in range(i, min(i + 2, len(target_labels))):
                                        target_labels[j] = formatted_del
                                        has_deletion = True
                                        processed_positions.add(del_start)
                                    break
                                else:
                                    target_labels[i] = formatted_del
                                    has_deletion = True
                                    processed_positions.add(del_start)
                                    break
                            elif del_count <= cut_pos - start:  # Add this condition for edited
                                target_labels[i] = formatted_del
                                is_edited = True
                                processed_positions.add(del_start)
                                break

                    for formatted_mis in formatted_mismatches:  # checking mismatches within targets
                        mis_count, mis_start = map(int, formatted_mis.split("I+"))
                        if start < mis_start < end_pos and (not has_deletion or has_deletion and mis_start not in processed_positions):
                            if target_labels[i] == "NONE":
                                target_labels[i] = formatted_mis
                                processed_positions.add(mis_start)
                            elif mis_start not in processed_positions:
                                target_labels[i] += f", {formatted_mis}"
                                processed_positions.add(mis_start)

                cellbarcode = read1_record.seq[:16]
                umi = read1_record.seq[16:]
                read1_header = read1_record.description.split()[0]
                row = [read2_header, cellbarcode, umi, ref_aligned, edited_aligned] + target_labels
                target_csv_writer.writerow(row)

                collapse_key = (cellbarcode, umi, edited_aligned)
                if collapse_key in collapsed_data:
                    collapsed_data[collapse_key][1].append(target_labels)
                else:
                    collapsed_data[collapse_key] = (read2_header, [target_labels])

        for (cellbarcode, umi, edited_aligned), (read2_header, target_labels_list) in collapsed_data.items():
            collapse_count = len(target_labels_list)
            common_target_labels = [max(targets) for targets in zip(*target_labels_list)]
            if collapse_count >= 1:
                collapse_row = [read2_header, cellbarcode, umi, collapse_count, ref_aligned, edited_aligned] + common_target_labels
                collapse_csv_writer.writerow(collapse_row)

           

if __name__ == "__main__":
    main()

