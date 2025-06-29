# scrna_barcode_tools_rs

A Rust utility to help grab cell-barcodes from FASTQ/BAM file. Current focused on PipSeq library kits, but can be extended for others.
This is not an official Illumina software utility, as is provided as-is 

# Assumptions:
- Illumina scRNA kit and use of DRAGEN to analyze
- Known Cell-barcodes file using this [format](https://help.dragen.illumina.com/product-guides/dragen-v4.3/dragen-single-cell-pipeline/dragen-scrna#known-barcode-lists)
- Cell barcode comes from the R1 or R2 FASTQ

# Installation instructions:
```bash
cd {INSTALL_DIR}
cargo build --release
```

## Common use-case : When a read does not have a valid barcode associated to it, how well does the extracted barcode match the expected barcodes list?

### GOALS:
- identify reads (primary alignments and pass QC) with no XB tag --- this indicates the read has no valid barcode
- for each of the reads above, revist initial fastq file to extract best barcode match
	- print out metrics of the barcode match

#### STEP1: extract reads that pass QC, and are not SECONDARY,DUPS, or SUPPLEMENTARY ALIGNMENTS - example command line(s) below
- samtools-1.21/samtools view -F 512 -F 256 -F 2048 sample1.scRNA.bam >  sample1.scRNA.final.sam
- samtools-1.21/samtools view -F 512 -F 256 -F 2048 sample2.scRNA.bam > sample2.scRNA.final.sam

#### STEP2: only grab reads with no XB tag - example command line(s) below
- grep -v XB sample1.scRNA.final.sam  > sample1.scRNA.final.missing_barcodes.sam 
- grep -v XB sample2.scRNA.final.sam > sample2.scRNA.final.missing_barcodes.sam

#### STEP3: obtain READ NAMES - example command line(s) below
- awk '{print $1}' sample1.scRNA.final.missing_barcodes.sam > sample1.scRNA.final.missing_barcodes.read_names.txt
- awk '{print $1}' sample2.scRNA.final.missing_barcodes.sam > sample2.scRNA.final.missing_barcodes.read_names.txt

#### STEP4: subset initial FASTQ using [seqtk](https://github.com/lh3/seqtk) - example command line(s) below
- seqtk subseq sample1_S2_L001_R1_001.fastq.gz sample1.scRNA.final.missing_barcodes.read_names.txt > sample1_S2_L001_R1_001.subsampled.fastq.gz

#### STEP5: for each read in the subsetted FASTQ
- extract barcode and indicate via hamming distance how close that barcode is to the whitelisted barcodes.
	- extracted barcodes and whitelisted barcodes are broken into blocks/tiers
	- for each block/tier we compute hamming distance and report back the hamming distance, extracted barcode block and best-matching whitelisted block

### command line to run
```bash
get_scrna_barcodes --fastq {FASTQ_path} --barcode-file {barcode_file} >  {OUTPUT_CSV}
```

![command line help](https://github.com/keng404/scrna_barcode_tools_rs/blob/main/Screen%20Shot%202025-06-25%20at%209.16.18%20AM.png)

### description of output csv ```{OUTPUT_CSV}``` file
- read_name --- read name that can be used to query FASTQ or BAM file
- barcode_sequence_extracted --- the entire barcode extracted from the read
- closest_whitelist_barcode --- the closest matching pipeseq barcode
- block1_dist --- hamming distance of the first block of the barcode compared to all Block1 of the pipseq barcodes
- block1_str --- barcode sequence of the first block of the barcode compared to all Block1 of the pipseq barcodes
- block2_dist --- hamming distance of the second block of the barcode compared to all Block2 of the pipseq barcodes
- block2_str --- barcode sequence of the second block of the barcode compared to all Block2 of the pipseq barcodes
- block3_dist --- hamming distance of the third block of the barcode compared to all Block3 of the pipseq barcodes
- block3_str --- barcode sequence of the third block of the barcode compared to all Block3 of the pipseq barcodes
- block4_dist --- hamming distance of the fourth block of the barcode compared to all Block4 of the pipseq barcodes
- block4_dist --- barcode sequence of the fourth block of the barcode compared to all Block4 of the pipseq barcodes
- total_dist --- total hamming distance of barcode_sequence_extracted and closest_whitelist_barcode

# DISCLAIMER

TO THE GREATEST EXTENT PERMITTED BY APPLICABLE LAW, THIS WEBSITE AND ITS CONTENT, INCLUDING ALL SOFTWARE, SOFTWARE CODE, SITE-RELATED SERVICES, AND DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. ALL WARRANTIES ARE REJECTED AND DISCLAIMED. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THE FOREGOING. ILLUMINA IS NOT OBLIGATED TO PROVIDE ANY SUPPORT FOR ANY OF THE FOREGOING, AND ANY SUPPORT ILLUMINA DOES PROVIDE IS SIMILARLY PROVIDED WITHOUT REPRESENTATION OR WARRANTY OF ANY KIND. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A REPRESENTATION OR WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY ILLUMINA.