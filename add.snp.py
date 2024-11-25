import gzip

# Step 1: Read the SNP file and create a dictionary
snp_dict = {}
with gzip.open('../hg19.snp.rsid.txt.gz', 'rt') as snp_file:
    for line in snp_file:
        parts = line.strip().split()
        chr_pos = (parts[0], parts[1])
        snp_dict[chr_pos] = parts[2]

# Function to process each file
def add_snp_column(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        header = infile.readline().strip() + '\tSNP\n'
        outfile.write(header)
        for line in infile:
            parts = line.strip().split()
            chr_pos = (f'chr{parts[0]}', parts[1])
            snp = snp_dict.get(chr_pos, 'NA')
            if snp != 'NA':
                outfile.write(line.strip() + f'\t{snp}\n')

# Step 2: Process each of the six files
files = ['BMI.UKBB.QC', 'BMI.BBJ.QC', 'TG.UKBB.QC', 'TG.BBJ.QC', 'Height.UKBB.QC', 'Height.BBJ.QC']
for file in files:
    add_snp_column(file, f'{file}_with_SNP.txt')
