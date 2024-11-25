import pandas as pd
import numpy as np
from scipy import stats

# Define sample sizes and corresponding file prefixes
sample_sizes = [
    {"bbj": 159095, "ukbb": 360388, "trait": "Height"},
    {"bbj": 105597, "ukbb": 343992, "trait": "TG"},
    {"bbj": 158284, "ukbb": 359983, "trait": "BMI"}
]

# Function to add Z_adj, Z_adj_norm, and P_adj columns
def add_adjusted_columns(df, sample_size):
    df['Z'] = pd.to_numeric(df['Z'], errors='coerce')  # Convert Z to numeric, set errors to NaN
    df = df.dropna(subset=['Z'])  # Drop rows where Z is NaN
    df['Z_adj'] = df['Z'] * np.sqrt(sample_size)
    df['Z_adj_norm'] = (df['Z_adj'] - df['Z_adj'].mean()) / df['Z_adj'].std()
    df['P_adj'] = 2 * (1 - stats.norm.cdf(np.abs(df['Z_adj_norm'])))
    return df

# Function to read a file with either space or tab separator
def read_file(file_path):
    try:
        return pd.read_csv(file_path, sep='\t', dtype=str)
    except pd.errors.ParserError:
        return pd.read_csv(file_path, sep=r'\s+', dtype=str)

# Process each combination of sample sizes and file prefixes
for sizes in sample_sizes:
    bbj_sample_size = sizes["bbj"]
    ukbb_sample_size = sizes["ukbb"]
    trait = sizes["trait"]
    
    # Load the data from the files for UKBB
    try:
        qc_with_snp_df_ukbb = read_file(f'{trait}.UKBB.QC_with_SNP.txt')
        
        # Add adjusted columns for UKBB sample size
        ukbb_df = add_adjusted_columns(qc_with_snp_df_ukbb.copy(), ukbb_sample_size)
        ukbb_df.to_csv(f'{trait}.UKBB.QC_with_SNP_{ukbb_sample_size}.txt', sep='\t', index=False)
        
        # Load the clumped file for UKBB
        clumped_df_ukbb = read_file(f'{trait}.UKBB.clumped')
        
        # Ensure the clumped file has at least three columns
        if clumped_df_ukbb.shape[1] >= 3:
            # Extract the third column and give it the header 'SNP'
            clumped_df_ukbb = clumped_df_ukbb.iloc[:, [2]]
            clumped_df_ukbb.columns = ['SNP']
            
            # Extract valid SNPs from the clumped file for UKBB
            valid_snps_ukbb = clumped_df_ukbb['SNP'].tolist()
            
            # Filter rows in QC_with_SNP where SNP is in valid_snps for UKBB
            filtered_qc_with_snp_df_ukbb = ukbb_df[ukbb_df['SNP'].isin(valid_snps_ukbb)]
            
            # Save the filtered dataframe
            filtered_qc_with_snp_df_ukbb.to_csv(f'{trait}.UKBB.filtered_{ukbb_sample_size}.txt', sep='\t', index=False)
        else:
            print(f"Clumped file for {trait} UKBB does not have enough columns.")
    
    except FileNotFoundError as e:
        print(f"File not found: {e.filename}")
    
    # Load the data from the files for BBJ
    try:
        qc_with_snp_df_bbj = read_file(f'{trait}.BBJ.QC_with_SNP.txt')
        
        # Add adjusted columns for BBJ sample size
        bbj_df = add_adjusted_columns(qc_with_snp_df_bbj.copy(), bbj_sample_size)
        bbj_df.to_csv(f'{trait}.BBJ.QC_with_SNP_{bbj_sample_size}.txt', sep='\t', index=False)
        
        # Load the clumped file for BBJ
        clumped_df_bbj = read_file(f'{trait}.BBJ.clumped')
        
        # Ensure the clumped file has at least three columns
        if clumped_df_bbj.shape[1] >= 3:
            # Extract the third column and give it the header 'SNP'
            clumped_df_bbj = clumped_df_bbj.iloc[:, [2]]
            clumped_df_bbj.columns = ['SNP']
            
            # Extract valid SNPs from the clumped file for BBJ
            valid_snps_bbj = clumped_df_bbj['SNP'].tolist()
            
            # Filter rows in QC_with_SNP where SNP is in valid_snps for BBJ
            filtered_qc_with_snp_df_bbj = bbj_df[bbj_df['SNP'].isin(valid_snps_bbj)]
            
            # Save the filtered dataframe
            filtered_qc_with_snp_df_bbj.to_csv(f'{trait}.BBJ.filtered_{bbj_sample_size}.txt', sep='\t', index=False)
        else:
            print(f"Clumped file for {trait} BBJ does not have enough columns.")
    
    except FileNotFoundError as e:
        print(f"File not found: {e.filename}")
