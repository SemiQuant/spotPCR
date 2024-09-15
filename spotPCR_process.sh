#!/bin/bash

# Command line inputs
while [[ $# -gt 0 ]]; do
  case $1 in
    -R1) R1="$2"; shift 2 ;;
    -R2) R2="$2"; shift 2 ;;
    -Ref_name) Ref_name="$2"; shift 2 ;;
    -sample_name) sample_name="$2"; shift 2 ;;
    -threads) threads="$2"; shift 2 ;;
    -umi) umi="$2"; shift 2 ;;
    -fqc) fqc="T"; shift 2 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

# Set defaults
bc="(?P<umi_1>.{10})T{2}.*"
# bc="(?P<umi_1>[AG]{1}.[CT]{1}.[ATCG]{3}.[AG]{1}.[CT]{1}.[ATCG]{3}.T{2}.*)"
umi="${umi:-$bc}"

# Check if required files exist
for file in "$R1" "$R2" "$Ref_name"; do
  if [ ! -f "$file" ]; then
    echo "Error: Required file $file is missing."
    exit 1
  fi
done

# Set default value for threads if not specified
threads="${threads:-4}"

# Default sample name if not specified
sample_name="${sample_name:-${R1/_R1*/}}"

# Run tools
if [ $fqc == "T" ]
then
  fastqc -q -t "$threads" "$R1"
  fastqc -q -t "$threads" "$R2"
fi

# these dont rerun unle
# Check if samtools index exists, if not, run samtools faidx
if [ ! -f "${Ref_name}.fai" ]
then
  samtools faidx "$Ref_name"
fi

# Check if BWA index files exist, if not, run bwa index
if [ ! -f "${Ref_name}.bwt" ] || [ ! -f "${Ref_name}.pac" ] || [ ! -f "${Ref_name}.ann" ] || [ ! -f "${Ref_name}.amb" ] || [ ! -f "${Ref_name}.sa" ]
then
  bwa index "$Ref_name"
fi

umi_tools extract \
  --extract-method regex \
  --bc-pattern "" \
  --bc-pattern2 "$umi" \
  -I "$R1" \
  --read2-in "$R2" \
  -S "${R1/.f*/.filt.fq}" \
  --read2-out "${R2/.f*/.filt.fq}" \
  -L "${sample_name}_extract.log"

bwa mem -t "$threads" "$Ref_name" "${R1/.f*/.filt.fq}" "${R2/.f*/.filt.fq}" > "${sample_name}.sam"
samtools sort -@ "$threads" "${sample_name}.sam" -o "${sample_name}_sorted.bam"
samtools index "${sample_name}_sorted.bam"
samtools flagstat "${sample_name}_sorted.bam" > "${sample_name}.flagstat.txt"


# umi_tools group --paired --per-gene --per-contig -I "${sample_name}_sorted.bam" --group-out="${sample_name}_grouped.tsv"
# # 5 or more umi per read
# awk -F'\t' '$6 >= 5 {print $4, $5, $6}' "${sample_name}_grouped.tsv" | sort -u | awk '{sum[$1] += $3} END {for (gene in sum) print gene, sum[gene]}' | sort -k2,2nr > "${sample_name}_usable_umis_reads.tsv"
# samtools idxstats "${sample_name}_sorted.bam" | awk '{OFS="\t"; if($1 != "*") print $1, $3}' > "${sample_name}_total_counts.txt"
# unique_umis=$(awk -F'\t' '$6 >= 5 {print $5}' "${sample_name}_grouped.tsv" | sort | uniq | wc -l)


# awk -v unique_umis="$unique_umis" '
# BEGIN {
#     OFS="\t"
#     print "Gene", "Usable_UMIs_Reads", "Total_Counts", "Unique_UMIs", "Percentage"
# }

# # Process the first file (total counts)
# NR==FNR {
#     counts[$1] = $2
#     next
# }

# # Process the second file (usable UMIs reads) and combine with total counts
# {
#     if ($1 in counts) {
#         percentage = (counts[$1] > 0) ? sprintf("%.2f", ($2 / counts[$1] * 100)) : "0.00"
#         print $1, $2, counts[$1], unique_umis, percentage
#     }
#     else {
#         print $1, $2, "0", unique_umis, "N/A"
#     }
# }' "${sample_name}_total_counts.txt" "${sample_name}_usable_umis_reads.tsv" > "${sample_name}_usable_umis_stats.tsv"



umi_tools dedup \
  --stdin="${sample_name}_sorted.bam" \
  --paired \
  --log="${sample_name}_dedup.log" \
  --output-stats="${sample_name}_dedup" > "${sample_name}_deduplicated_reads.bam"

samtools index "${sample_name}_deduplicated_reads.bam"

freebayes -p 1 -f "$Ref_name" "${sample_name}_sorted.bam" > "${sample_name}_orig.vcf"
freebayes -p 1 -f "$Ref_name" "${sample_name}_deduplicated_reads.bam" > "${sample_name}.vcf"

vcf_file="${sample_name}.vcf"
vcf_orig_file="${sample_name}_orig.vcf"

python3 -c '
import vcf
import pandas as pd
import os

def vcf_to_tsv(vcf_input):
    output_file = os.path.splitext(vcf_input)[0] + ".tsv"
    vcf_reader = vcf.Reader(open(vcf_input, "r"))
    df = pd.DataFrame(
        [
            (
                record.POS,
                record.CHROM,
                record.REF,
                record.ALT[0] if record.ALT else "",
                record.INFO.get("DP", 0),
                record.INFO.get("AO", [0])[0] if record.INFO.get("AO") else 0,
                record.INFO.get("AC", [0])[0] if record.INFO.get("AC") else 0,
                round(record.INFO.get("AO", [0])[0] / record.INFO.get("DP", 0), 2) if record.INFO.get("DP") != 0 else 0,
                record.QUAL if record.QUAL else 0
            )
            for record in vcf_reader
        ],
        columns=["POS", "CHROM", "REF", "ALT", "Depth", "Alt Depth", "AC", "Proportion Depth", "QUAL"]
    )
    df.to_csv(output_file, sep="\t", index=False)
    print(f"Saved {vcf_input} as {output_file}")

# Process both VCF files
vcf_to_tsv("'"$vcf_file"'")
vcf_to_tsv("'"$vcf_orig_file"'")
'

bam="${sample_name}_deduplicated_reads.bam"
REGIONS=$(samtools view -H "$bam" | grep "@SQ" | cut -f 2 | cut -d ':' -f 2)

for REGION in $REGIONS; do
    for BAM_FILE in "${sample_name}_sorted.bam" "${sample_name}_deduplicated_reads.bam"; do
        DEPTH_VALUES=$(samtools depth -a -r "$REGION" "$BAM_FILE" | cut -f3)
        AVG_DEPTH=$(echo "$DEPTH_VALUES" | awk '{ total += $1; count++ } END { print total/count }')
        echo "Bam file $BAM_FILE Region $REGION: Average Depth = $AVG_DEPTH" >> "${sample_name}.stats.txt"
        umi_tools count --per-gene --gene-tag=pncA --per-cell -I $BAM_FILE -S "${sample_name}_${REGION}.counts.tsv"
    done
done

###
# # Calculate total reads per contig
# total_reads=$(samtools idxstats "${sample_name}_sorted.bam" | awk '{OFS="\t"; if($1 != "*") sum+=$3} END {print sum}')

# # Calculate usable reads (UMIs observed 5 or more times)
# usable_reads=$(awk 'NR>1 && $2 >= 5 {sum+=$4} END {print sum}' "${sample_name}_dedup_per_umi.tsv")

# # Calculate unique UMIs (observed 5 or more times)
# unique_umis=$(awk 'NR>1 && $2 >= 5 {count++} END {print count}' "${sample_name}_dedup_per_umi.tsv")

# echo -e "total_reads\tusable_reads\tunique_umis" > "${sample_name}_basicstats.tsv"
# echo -e "${total_reads}\t${usable_reads}\t${unique_umis}" >> "${sample_name}_basicstats.tsv"
###

# bam="${sample_name}_deduplicated_reads.bam"

# umi_tools count --extract-umi-method=read_id --umi-separator="_" --per-contig --per-gene --wide-format-cell-counts -I "$bam" -S "${sample_name}.counts.tsv"
# samtools idxstats "$bam" | awk '{OFS="\t"; if($1 != "*") print $1, $3}' > "${sample_name}.total_counts.txt"

# # Combine UMI counts and total counts
# echo -e "contig\tumi_count\ttotal_reads" > "${sample_name}combined_counts.tsv"
# join -t $'\t' <(sort "${sample_name}.counts.tsv") <(sort "${sample_name}.total_counts.txt") >> "${sample_name}combined_counts.tsv"


# Cleanup
rm "${sample_name}.sam" "${R2/.f*/.filt.fq}" "${R1/.f*/.filt.fq}" 
#"${sample_name}_usable_umis_reads.tsv" "${sample_name}_total_counts.txt"
# "${sample_name}.counts.tsv" "${sample_name}.total_counts.txt"

exit 0
