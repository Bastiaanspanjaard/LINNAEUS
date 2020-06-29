#!/bin/bash
# Count reads per barcode and select barcodes with more than 500 reads as cellular barcodes
awk 'NR % 4 == 2 {printf("%s\n", substr($1, 1, 16))}' P5_scar_R1.fastq | sort | uniq -c > P5_scar_bc_freqs.txt
awk 'BEGIN{OFS = "\t"; bc = 1} $1 > 500 {print bc, $2; bc++}' P5_scar_bc_freqs.txt > P5_scar_barcodes.csv
# Create UMIBC-file
awk '(NR % 4 == 1) {print $1, $2} (NR % 4 == 3) {print $1} (NR % 2 == 0) {printf("%s%s\n", substr($1, 17, 10), substr($1, 1, 16))}' P5_scar_R1.fastq > P5_scar_UMIBC.fastq
# Map and extract scars
/local/Bastiaan/Scripts/scar_CIGAR_sc_10X_v2.pl -R1=P5_scar_UMIBC.fastq -R2=P5_scar_R2.fastq -op=P5_scar -t=1 -r=/local/gene_models/lintrace_hGFP_RFP_ERCC92.fa -bc=P5_scar_barcodes.csv -g=RFP -k=75 -mU=1 -l=75 -ps=GAGTTCAAGACCATCTACATGGCC
# Count how often every scar is sequenced, remove those that only occur once, filter the rest
tail -n +2 P5_scar_scars.txt | sort -k2,2 -k3,3 -k1,1 -k7,7 | uniq -c | cut -f1-4,7 > P5_scar_reads.txt
awk '$1>1' P5_scar_reads.txt > P5_scar_reads_over1.txt
source /local/Bastiaan/Scripts/scar_filter/bin/activate
/local/Bastiaan/bitbucket_scripts/sc_scar/scar_filter.py -i P5_scar_reads_over1.txt -o P5_scar_filtered_scars.csv
