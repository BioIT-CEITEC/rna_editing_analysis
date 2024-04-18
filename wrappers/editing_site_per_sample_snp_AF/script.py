######################################
# wrapper for rule: per_sample_snp_AF_computing
######################################
from snakemake.shell import shell
import csv

log_filename = str(snakemake.log)

def save_sorted_table(file_path, output_file_path):
    # Load data from the file
    data = []
    with open(file_path, newline='') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            if len(row) >= 2:  # Check if at least two columns are present
                data.append((row[0], int(row[1])))

    # Sort data by chromosome and position
    data.sort(key=lambda x: (x[0], x[1]))

    # Write the sorted data to a new file
    with open(output_file_path, 'w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        writer.writerow(["Chromosome", "Position"])  # Optional: write a header
        for chrom, pos in data:
            writer.writerow([chrom, pos])

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: per_sample_snp_AF_computing \n##\n")
f.close()

# File path for the sorted output
edit_pos_filename = snakemake.input.snp_tsv

# Call the function to save the sorted table
save_sorted_table(snakemake.input.snp_tsv, edit_pos_filename.replace(".tsv", ".alleleCounterFormat.tsv"))

command = "alleleCounter --dense-snps" + \
            " -r " + snakemake.input.ref + \
            " -l " + edit_pos_filename.replace(".tsv", ".alleleCounterFormat.tsv") + \
            " -b  " + snakemake.input.bam + \
            " -o " + snakemake.output.snp_tab + \
            " >> " + log_filename + " 2>&1 "

f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()

shell(command)