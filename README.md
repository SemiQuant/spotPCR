This Bash script performs UMI identification, read alignment, deduplication, and variant calling for spotPCR sequence data. 

**Usage**

To use this script, execute it from the command line with the following command:

```
./spotPCR_process.sh \
  -R1 <Read1.fastq> \
  -R2 <Read2.fastq> \
  -Ref_name <Reference.fasta> \
  -sample_name <SampleName> \
  -threads <Threads> \
  -umi <UMI_Barcode_Pattern>
```

**Command-line Arguments**

| Flag         | Description                                                                                                            |
|--------------|------------------------------------------------------------------------------------------------------------------------|
| -R1          | Path to the first input read file (Read 1).                                                                            |
| -R2          | Path to the second input read file (Read 2).                                                                           |
| -Ref_name    | Path to the reference genome or sequence file.                                                                         |
| -sample_name | Name to be used for the output files (optional).                                                                       |
| -threads     | Number of threads for parallel processing (optional, default is 4).                                                    |
| -umi         | UMI (Unique Molecular Identifier) regex barcode pattern for read extraction (optional, default is a '(?P.{10})T{2}.*') |


**Requirements**

- [fastqc: Quality control tool for sequencing data.](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [samtools: Utilities for manipulating sequence data in SAM and BAM formats.](https://www.htslib.org/)
- [bwa: Burrows-Wheeler Aligner for read alignment.](https://bio-bwa.sourceforge.net/)
- [umi_tools: A tool for working with UMIs.](https://umi-tools.readthedocs.io/en/latest/QUICK_START.html)
- [freebayes: A Bayesian genetic variant detector.](https://github.com/freebayes/freebayes)
- python3: Python for data manipulation.
    - pyVCF
    - pandas


**Output**

The script generates various output files, including aligned BAM files, a VCF file containing variant calls, and a statistics file for depth analysis.


**Run Example**

./spotPCR_process.sh \
  -R1 ./example/ex_R1_001.fastq.gz \
  -R2 ./example/ex_R2_001.fastq.gz \
  -Ref_name ./reference/ref.fasta \
  -sample_name example \
  -threads 4
```

The expected outputs are in the example folder for comparison.