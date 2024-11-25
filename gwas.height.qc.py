import pandas as pd
import numpy as np

# Function to normalize Z-score
def normalize_z(z):
    return (z - np.mean(z)) / np.std(z)

# Process BBJ Dataset
def process_bbj(input_file, output_file):
    bbj_df = pd.read_csv(input_file, sep=' ', compression='gzip')
    bbj_df['Z'] = bbj_df['BETA'] / bbj_df['SE']
    bbj_df['Z_norm'] = normalize_z(bbj_df['Z'])
    bbj_df.rename(columns={'ALT_freq': 'Frq'}, inplace=True)
    bbj_df.rename(columns={'P_BOLT': 'P'}, inplace=True)
    bbj_df[['CHR', 'POS', 'ALT', 'Frq', 'Z', 'Z_norm', 'P']].to_csv(output_file, sep='\t', index=False)

# Process UKBB Dataset
def process_ukbb(input_file, output_file):
    ukbb_df = pd.read_csv(input_file, sep='\t', compression='gzip')
    ukbb_df[['CHR', 'POS', 'REF', 'ALT']] = ukbb_df['variant'].str.split(':', expand=True)
    ukbb_df['Frq'] = np.where(ukbb_df['ALT'] == ukbb_df['minor_allele'], ukbb_df['minor_AF'], 1 - ukbb_df['minor_AF'])
    ukbb_df['Z'] = ukbb_df['beta'] / ukbb_df['se']
    ukbb_df['Z_norm'] = normalize_z(ukbb_df['Z'])
    ukbb_df.rename(columns={'pval': 'P'}, inplace=True)
    ukbb_df[['CHR', 'POS', 'ALT', 'Frq', 'Z', 'Z_norm', 'P']].to_csv(output_file, sep='\t', index=False)

# File paths
bbj_input = './BBJ/height_autosomes_BOLT.txt.gz'
bbj_output = './QC.GWAS/Height.BBJ.QC'
ukbb_input = './UKBB/height_raw.gwas.imputed_v3.both_sexes.tsv.bgz'
ukbb_output = './QC.GWAS/Height.UKBB.QC'

# Run the processes
process_bbj(bbj_input, bbj_output)
process_ukbb(ukbb_input, ukbb_output)
