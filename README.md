# scrna_barcode_tools_rs

## steps taken to get to the output
# GOALS:
	- identify reads (primary alignments and pass QC) with no XB tag --- this indicates the read has no valid barcode
	- for each of the reads above, revist initial fastq file to extract best barcode match
		- print out metrics of the barcode match

### STEP1: extract reads that pass QC, and are not SECONDARY,DUPS, or SUPPLEMENTARY ALIGNMENTS
###samtools-1.21/samtools view -F 512 -F 256 -F 2048 XRQ-DM-032725.scRNA.bam >  XRQ-DM-032725.scRNA.final.sam
####samtools-1.21/samtools view -F 512 -F 256 -F 2048 PSC8-DM-032725.scRNA.bam > PSC8-DM-032725.scRNA.final.sam
#### STEP2: only grab reads with no XB tag
##### grep -v XB XRQ-DM-032725.scRNA.final.sam  > XRQ-DM-032725.scRNA.final.missing_barcodes.sam 
######  grep -v XB PSC8-DM-032725.scRNA.final.sam > PSC8-DM-032725.scRNA.final.missing_barcodes.sam
###### STEP3: obtain READ NAMES
#### awk '{print $1}' PSC8-DM-032725.scRNA.final.missing_barcodes.sam > $HOME/PSC8-DM-032725.scRNA.final.missing_barcodes.read_names.txt
### awk '{print $1}' XRQ-DM-032725.scRNA.final.missing_barcodes.sam > $HOME/XRQ-DM-03272.scRNA.final.missing_barcodes.read_names.txt
#####
##### STEP4: subset initial FASTQ using seqtk
### seqtk subseq /data/PSC8-DM-032725_ds.41dfc6fa977842cc98ef01f66500d6be/PSC8-DM-032725_S2_L001_R1_001.fastq.gz /additional_data/PSC8-DM-032725.scRNA.final.missing_barcodes.read_names.txt > /data/PSC8-DM-032725_S2_L001_R1_001.subsampled.fastq.gz
####### STEP5: for each read in the subsetted FASTQ
extract barcode and indicate via hamming distance how close that barcode is to the whitelisted barcodes.
	- extracted barcodes and whitelisted barcodes are broken into blocks/tiers
	- for each block/tier we compute hamming distance and report back the hamming distance, extracted barcode block and best-matching whitelisted block

### description of output csv file
read_name --- read name that can be used to query FASTQ or BAM file
barcode_sequence_extracted --- the entire barcode extracted from the read
closest_whitelist_barcode --- the closest matching pipeseq barcode
block1_dist --- hamming distance of the first block of the barcode compared to all Block1 of the pipseq barcodes
block1_str --- barcode sequence of the first block of the barcode compared to all Block1 of the pipseq barcodes
block2_dist --- hamming distance of the second block of the barcode compared to all Block2 of the pipseq barcodes
block2_str --- barcode sequence of the second block of the barcode compared to all Block2 of the pipseq barcodes
block3_dist --- hamming distance of the third block of the barcode compared to all Block3 of the pipseq barcodes
block3_str --- barcode sequence of the third block of the barcode compared to all Block3 of the pipseq barcodes
block4_dist --- hamming distance of the fourth block of the barcode compared to all Block4 of the pipseq barcodes
block4_dist --- barcode sequence of the fourth block of the barcode compared to all Block4 of the pipseq barcodes
total_dist --- total hamming distance of barcode_sequence_extracted and closest_whitelist_barcode
