import pandas as pd
import numpy as np
import scipy.stats as stats

# Function to calculate Z-score
def calculate_z(beta, se):
    return beta / se

# Function to adjust Z-score by sample size
def adjust_z(z, n):
    return z * np.sqrt(n)

# Function to calculate p-value from Z-score
def calculate_p_value(z):
    return 2 * (1 - stats.norm.cdf(abs(z)))

# Process BBJ Dataset
def process_bbj(input_file, output_file, sample_size):
    bbj_df = pd.read_csv(input_file, sep='\t', compression='gzip')
    bbj_df['Z'] = calculate_z(bbj_df['BETA'], bbj_df['SE'])
    bbj_df['Z_adjusted'] = adjust_z(bbj_df['Z'], sample_size)
    bbj_df['P_adjusted'] = bbj_df['Z_adjusted'].apply(calculate_p_value)
    bbj_df['P_adjusted'] = bbj_df['P_adjusted'].apply(lambda x: f"{x:.2e}")
    bbj_df[['CHR', 'POS', 'ALT', 'Frq', 'Z', 'Z_adjusted', 'P_adjusted']].to_csv(output_file, sep='\t', index=False)

# Process UKBB Dataset
def process_ukbb(input_file, output_file, sample_size):
    ukbb_df = pd.read_csv(input_file, sep='\t', compression='gzip')
    ukbb_df[['CHR', 'POS', 'REF', 'ALT']] = ukbb_df['variant'].str.split(':', expand=True)
    ukbb_df['Frq'] = np.where(ukbb_df['ALT'] == ukbb_df['minor_allele'], ukbb_df['minor_AF'], 1 - ukbb_df['minor_AF'])
    ukbb_df['Z'] = calculate_z(ukbb_df['beta'], ukbb_df['se'])
    ukbb_df['Z_adjusted'] = adjust_z(ukbb_df['Z'], sample_size)
    ukbb_df['P_adjusted'] = ukbb_df['Z_adjusted'].apply(calculate_p_value)
    ukbb_df['P_adjusted'] = ukbb_df['P_adjusted'].apply(lambda x: f"{x:.2e}")
    ukbb_df.rename(columns={'pval': 'P'}, inplace=True)
    ukbb_df[['CHR', 'POS', 'ALT', 'Frq', 'Z', 'Z_adjusted', 'P_adjusted']].to_csv(output_file, sep='\t', index=False)

# File paths and sample sizes
bbj_input = './BBJ/TG.autosome.txt.gz'
bbj_output = './adj.GWAS/TG.BBJ.adj'
ukbb_input = './UKBB/TG_raw.gwas.imputed_v3.both_sexes.varorder.tsv.bgz'
ukbb_output = './adj.GWAS/TG.UKBB.adj'
bbj_sample_size = 105597  # Example sample size for BBJ
ukbb_sample_size = 343992  # Example sample size for UKBB

# Run the processes
process_bbj(bbj_input, bbj_output, bbj_sample_size)
process_ukbb(ukbb_input, ukbb_output, ukbb_sample_size)
