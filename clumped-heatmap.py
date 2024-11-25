import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Function to process and match datasets
def process_and_match(bbj_file, ukbb_file, output_file):
    # Read BBJ and UKBB datasets
    bbj_df = pd.read_csv(bbj_file, sep='\t')
    ukbb_df = pd.read_csv(ukbb_file, sep='\t')
    
    # Ensure the 'CHR' column is of the same type in both dataframes
    bbj_df['CHR'] = bbj_df['CHR'].astype(str)
    ukbb_df['CHR'] = ukbb_df['CHR'].astype(str)
    
    # Merge datasets on CHR and POS
    merged_df = pd.merge(bbj_df, ukbb_df, on=['CHR', 'POS'], suffixes=('_BBJ', '_UKBB'))
    
    # Check conditions and adjust UKBB columns
    condition = merged_df['ALT_BBJ'] != merged_df['ALT_UKBB']
    merged_df.loc[condition, 'Frq_UKBB'] = 1 - merged_df['Frq_UKBB']
    merged_df.loc[condition, 'Z_norm_UKBB'] = -merged_df['Z_norm_UKBB']
    
    # Save the merged data without filtering
    merged_df.to_csv(output_file, sep='\t', index=False)
    
    return merged_df

# Function to plot heatmap with enlarged axis label font size and specified titles
def plot_heatmap(df, bbj_title, ukbb_title, output_file):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6), sharey=True)
    
    # Determine the color scale limits
    vmin = 0  # Change the color interval to >0
    vmax = max(df['Z_norm_BBJ'].max(), df['Z_norm_UKBB'].max())
    
    # Plot BBJ heatmap with legend
    sns.histplot(data=df, x='Frq_BBJ', y='Z_norm_BBJ', bins=30, cbar=True, color='orange', ax=ax1, vmin=vmin, vmax=vmax)
    ax1.set_xlim(0, 1)  # Adjust x-axis limits to match UKBB heatmap
    ax1.set_title(bbj_title)
    ax1.set_xlabel('Effect allele frequency', fontsize=14)
    ax1.set_ylabel('Normalized sample size weighted Z score', fontsize=14)
    
    # Plot UKBB heatmap with legend limited to values >0
    sns.histplot(data=df, x='Frq_UKBB', y='Z_norm_UKBB', bins=30, cbar=True, color='orange', ax=ax2, vmin=vmin, vmax=vmax)
    ax2.set_xlim(0, 1)  # Adjust x-axis limits to match BBJ heatmap
    ax2.set_title(ukbb_title)
    ax2.set_xlabel('Effect allele frequency', fontsize=14)
    ax2.set_ylabel('')
    
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

# File pairs
file_pairs = [
    ('TG.BBJ.QC.5e8_105597.txt', 'TG.UKBB.QC.5e8_343992.txt', 'TG.matched.QC.txt', 'TG.BBJ', 'TG.UKBB'),
    ('Height.BBJ.QC.5e8_159095.txt', 'Height.UKBB.QC.5e8_360388.txt', 'Height.matched.QC.txt', 'Height.BBJ', 'Height.UKBB'),
    ('BMI.BBJ.QC.5e8_158284.txt', 'BMI.UKBB.QC.5e8_359983.txt', 'BMI.matched.QC.txt', 'BMI.BBJ', 'BMI.UKBB')
]

# Process each pair and create heatmaps
for bbj_file, ukbb_file, matched_file, bbj_title, ukbb_title in file_pairs:
    matched_df = process_and_match(bbj_file, ukbb_file, matched_file)
    plot_heatmap(matched_df, bbj_title, ukbb_title, f'{bbj_title}.heatmap.pdf')

print("Processing complete. Merged files and heatmaps are saved.")
