import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Function to normalize Z-score
def normalize_z(z):
    return (z - np.mean(z)) / np.std(z)

# Function to process and match datasets
def process_and_match(bbj_file, ukbb_file, output_file):
    # Read BBJ and UKBB datasets
    bbj_df = pd.read_csv(bbj_file, sep='\t')
    ukbb_df = pd.read_csv(ukbb_file, sep='\t')
    
    # Merge datasets on CHR and POS
    merged_df = pd.merge(bbj_df, ukbb_df, on=['CHR', 'POS'], suffixes=('_BBJ', '_UKBB'))
    
    # Check conditions and adjust UKBB columns
    condition = merged_df['ALT_BBJ'] != merged_df['ALT_UKBB']
    merged_df.loc[condition, 'Frq_UKBB'] = 1 - merged_df['Frq_UKBB']
    merged_df.loc[condition, 'Z_norm_UKBB'] = -merged_df['Z_norm_UKBB']
    
    # Filter based on P values
    filtered_df = merged_df[(merged_df['P_BBJ'] < 5e-3) & (merged_df['P_UKBB'] < 5e-3)]
    
    # Save the filtered data
    filtered_df.to_csv(output_file, sep='\t', index=False)
    
    return filtered_df

# Function to plot heatmap
def plot_heatmap(df, bbj_title, ukbb_title, output_file):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6), sharey=True)
    
    # Determine the color scale limits
    vmin = min(df['Z_norm_BBJ'].min(), df['Z_norm_UKBB'].min())
    vmax = max(df['Z_norm_BBJ'].max(), df['Z_norm_UKBB'].max())
    
    # Plot BBJ heatmap
    sns.histplot(data=df, x='Frq_BBJ', y='Z_norm_BBJ', bins=30, cbar=True, color='orange', ax=ax1, vmin=vmin, vmax=vmax)
    ax1.set_title(bbj_title)
    ax1.set_xlabel('Frq')
    ax1.set_ylabel('Z_norm')
    
    # Plot UKBB heatmap
    sns.histplot(data=df, x='Frq_UKBB', y='Z_norm_UKBB', bins=30, cbar=True, color='orange', ax=ax2, vmin=vmin, vmax=vmax)
    ax2.set_title(ukbb_title)
    ax2.set_xlabel('Frq')
    ax2.set_ylabel('')
    
    # Add a single legend
    plt.legend(['BBJ', 'UKBB'], loc='upper right')
    
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

# File pairs
file_pairs = [
    ('./QC.GWAS/Height.BBJ.QC', './QC.GWAS/Height.UKBB.QC', './QC.GWAS/Height.matched.QC', 'Height'),
    ('./QC.GWAS/TG.BBJ.QC', './QC.GWAS/TG.UKBB.QC', './QC.GWAS/TG.matched.QC', 'TG'),
    ('./QC.GWAS/HDL.BBJ.QC', './QC.GWAS/HDL.UKBB.QC', './QC.GWAS/HDL.matched.QC', 'HDL'),
    ('./QC.GWAS/BMI.BBJ.QC', './QC.GWAS/BMI.UKBB.QC', './QC.GWAS/BMI.matched.QC', 'BMI')
]

# Process each pair and create heatmaps
for bbj_file, ukbb_file, matched_file, title in file_pairs:
    matched_df = process_and_match(bbj_file, ukbb_file, matched_file)
    plot_heatmap(matched_df, f'{title}.BBJ', f'{title}.UKBB', f'./QC.GWAS/{title}.heatmap.pdf')
