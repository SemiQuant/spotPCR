#!/bin/python
import pysam
import collections
import re
import pandas as pd
import numpy as np
import os
import argparse

parser = argparse.ArgumentParser(description="Inputs")
parser.add_argument("-d", "--directory", 
                    default=os.getcwd(),
                    help="Path to the directory to process (default: current directory)")
parser.add_argument("-u", "--min_umi_count", 
                    type=int, 
                    default=50,
                    help="Minimum UMI count (default: 50)")
args = parser.parse_args()

bam_directory = args.directory
min_umi_count = args.min_umi_count

def process_bam_file(file_path):
    bam_file = pysam.AlignmentFile(file_path, "rb")
    umi_counts = collections.defaultdict(lambda: collections.defaultdict(int))
    read_counts = collections.defaultdict(int)
    unique_umis = collections.defaultdict(set)
    total_reads = aligned_reads = reads_with_umi = 0
    umi_pattern = re.compile(r'_([ACGT]{10})$')

    for read in bam_file:
        total_reads += 1
        if not read.is_unmapped:
            aligned_reads += 1
            umi_match = umi_pattern.search(read.query_name)
            if umi_match:
                reads_with_umi += 1
                ref_name = read.reference_name
                umi = umi_match.group(1)
                umi_counts[ref_name][umi] += 1
                unique_umis[ref_name].add(umi)

    total_usable_reads = sum(count for ref_umis in umi_counts.values() 
                             for count in ref_umis.values() if count >= min_umi_count)
    for ref_name, umis in umi_counts.items():
        read_counts[ref_name] = sum(count for count in umis.values() if count >= min_umi_count)

    data = [{'Filename': os.path.basename(file_path),
             'CHROM': ref_name,
             'Reads with usable UMIs': count,
             'Total aligned reads': aligned_reads,
             'Unique UMIs': len(unique_umis[ref_name]),
             'Proportion': count / aligned_reads if aligned_reads > 0 else 0}
            for ref_name, count in read_counts.items()]

    total_unique_umis = sum(len(umis) for umis in unique_umis.values())
    df = pd.DataFrame(data + [{
        'Filename': os.path.basename(file_path),
        'CHROM': 'Total',
        'Reads with usable UMIs': total_usable_reads,
        'Total aligned reads': aligned_reads,
        'Unique UMIs': total_unique_umis,
        'Proportion': total_usable_reads / aligned_reads if aligned_reads > 0 else 0
    }])

    bam_file.close()
    return df

os.chdir(bam_directory)
bam_files = [f for f in os.listdir() if f.endswith('_sorted.bam')]

all_data = [process_bam_file(bam_file) for bam_file in bam_files]
combined_df = pd.concat(all_data, ignore_index=True)

for col in ['Total aligned reads', 'Reads with usable UMIs', 'Unique UMIs']:
    combined_df[col] = combined_df[col].apply(lambda x: f"{x:,}")
combined_df['Proportion'] = combined_df['Proportion'].apply(lambda x: f"{x:.2%}")

combined_df.to_csv('alignment_umi_stats.csv', index=False)

# combined_df['grp'] = combined_df['Filename']
combined_df['grp'] = combined_df['Filename'].apply(lambda x: re.search(r'.*_([A-Za-z0-9]+)_S.*', x).group(1) if re.search(r'.*_([A-Za-z0-9]+)_S.*', x) else None)

for col in ['Reads with usable UMIs', 'Total aligned reads', 'Unique UMIs']:
    combined_df[col] = pd.to_numeric(combined_df[col].str.replace(',', ''), errors='coerce')

combined_df['prop'] = combined_df['Reads with usable UMIs'] / combined_df['Total aligned reads'] * 100

def calculate_stats(group):
    return pd.Series({
        'prop_mean': group['prop'].mean(),
        'prop_sd': group['prop'].std(),
        'prop_median': group['prop'].median(),
        'prop_iqr': f"{group['prop'].quantile(0.75):.2f}, {group['prop'].quantile(0.25):.2f}",
        'unique_umis_mean': group['Unique UMIs'].mean(),
        'unique_umis_sd': group['Unique UMIs'].std(),
        'unique_umis_median': group['Unique UMIs'].median(),
        'unique_umis_iqr': f"{group['Unique UMIs'].quantile(0.75):.0f}, {group['Unique UMIs'].quantile(0.25):.0f}",
        'total_unique_umis': group['Unique UMIs'].sum()
    })

result = combined_df.groupby('grp').apply(calculate_stats).reset_index()

# Format the results
for col in result.columns:
    if col not in ['grp', 'prop_iqr', 'unique_umis_iqr']:
        result[col] = result[col].apply(lambda x: f"{x:.2f}")

print(result.to_string(index=False))