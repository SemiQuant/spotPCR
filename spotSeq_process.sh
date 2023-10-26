#!/bin/bash
# Command line inputs
while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
    -R1)
      R1="$2"
      shift # past argument
      shift # past value
      ;;
    -R2)
      R2="$2"
      shift
      shift
      ;;
    -Ref_name)
      Ref_name="$2"
      shift
      shift
      ;;
    -sample_name)
      sample_name="$2"
      shift
      shift
      ;;
    -threads)
      threads="$2"
      shift
      shift
      ;;
    -umi)
      umi="$2"
      shift
      shift
      ;;
    *)
      # unknown option
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

# set defaults
bc="(?P<umi_1>.{10})T{2}.*"
umi="${umi-$bc}"
# "(?P<umi_1>.{10})T{2}.*"


# Check if required files exist
if [ ! -f "$R1" ] || [ ! -f "$R2" ] || [ ! -f "$Ref_name" ]; then
  echo "Error: One or more required files are missing."
  echo "Please make sure the following files exist:"
  echo "R1: $R1"
  echo "R2: $R2"
  echo "Ref_name: $Ref_name"
  exit 1
fi

# Check if threads is still empty (user didn't specify it)
if [ -z "$threads" ]; then
  threads=4  # Assign the default value
fi

if [ -z "$sample_name" ]; then
  sample_name="${R1/_R1*/}"
fi


# fastqc -q -t "$threads" "$R1"
# fastqc -q -t "$threads" "$R2"

samtools faidx "$Ref_name"
bwa index "$Ref_name"

# bwa mem -t "$threads" "$Ref_name" "${R1}" "${R2}" > "${sample_name}.sam"
# samtools sort -@ "$threads" "${sample_name}.sam" -o "${sample_name}_sorted.bam"
# samtools index "${sample_name}_sorted.bam"


umi_tools extract \
  --extract-method regex \
  --bc-pattern "" \
  --bc-pattern2 "(?P<umi_1>.{10})T{2}.*" \
  -I "$R1" \
  --read2-in "$R2" \
  -S "${R1/.f*/.filt.fq}" \
  --read2-out "${R2/.f*/.filt.fq}" \
  -L "${sample_name}_extract.log"

# or
# --bc-pattern2 "(?P<discard_1>TGCGCCGAGCGTGAACTCAG){s<=1}(?P<umi_1>.{10})T{2}.*" \

bwa mem -t "$threads" "$Ref_name" "${R1/.f*/.filt.fq}" "${R2/.f*/.filt.fq}" > "${sample_name}.sam"
samtools sort -@ "$threads" "${sample_name}.sam" -o "${sample_name}_sorted.bam"
samtools index "${sample_name}_sorted.bam"

umi_tools dedup \
  --stdin="${sample_name}_sorted.bam" \
  --paired \
  --log="${sample_name}__dedup.log" \
  --output-stats="${sample_name}_dedup" > "${sample_name}_deduplicated_reads.bam"

freebayes -p 1 -f "$Ref_name" "${sample_name}_sorted.bam" > "${sample_name}.vcf"

vcf_file="${sample_name}.vcf" && python3 -c 'import vcf, pandas as pd; vcf_reader = vcf.Reader(open("'"$vcf_file"'", "r")); df = pd.DataFrame([(record.POS, record.REF, record.ALT[0] if record.ALT else "", record.INFO.get("DP", 0), record.INFO.get("AO", [0])[0] if record.INFO.get("AO") else 0, record.INFO.get("AC", [0])[0] if record.INFO.get("AC") else 0, round(record.INFO.get("AO", [0])[0] / record.INFO.get("DP", 0), 2) if record.INFO.get("DP") != 0 else 0, record.QUAL if record.QUAL else 0) for record in vcf_reader], columns=["POS", "REF", "ALT", "Depth", "Alt Depth", "AC", "Proportion Depth", "QUAL"]); print(df)'


for REGION in "${REGIONS[@]}"
do
    bam="${sample_name}_sorted.bam"
    REGIONS=$(samtools view -H "$bam" | grep "@SQ" | cut -f 2 | cut -d ':' -f 2)
    DEPTH_VALUES=$(samtools depth -a -r $REGION $BAM_FILE | cut -f3)
    AVG_DEPTH=$(echo "$DEPTH_VALUES" | awk '{ total += $1; count++ } END { print total/count }')
    echo "Bam file region $REGION: Average Depth = $AVG_DEPTH" > "${sample_name}.stats.tsv"

    bam="${sample_name}_deduplicated_reads.bam"
    DEPTH_VALUES=$(samtools depth -a -r $REGION $BAM_FILE | cut -f3)
    AVG_DEPTH=$(echo "$DEPTH_VALUES" | awk '{ total += $1; count++ } END { print total/count }')
    echo "Bam file deduplicated Region $REGION: Average Depth = $AVG_DEPTH" >> "${sample_name}.stats.tsv"
done


# Cleanup
rm "${sample_name}.sam" "${R2/.f*/.filt.fq}" "${R1/.f*/.filt.fq}"


