import pandas as pd
import glob

# Function to read a file with either space or tab separator
def read_file(file_path):
    try:
        return pd.read_csv(file_path, sep='\t', dtype=str)
    except pd.errors.ParserError:
        return pd.read_csv(file_path, sep=r'\s+', dtype=str)

# Extract rows from *.QC_with_SNP_*.txt files where column P_adj < 5e-8 and save to *.QC.5e8 files
for file_path in glob.glob('*QC_with_SNP_*.txt'):
    df = read_file(file_path)
    df_filtered = df[pd.to_numeric(df['P_adj'], errors='coerce') < 5e-8]
    output_file_path = file_path.replace('QC_with_SNP_', 'QC.5e8_')
    df_filtered.to_csv(output_file_path, sep='\t', index=False)

print("Extraction complete.")
