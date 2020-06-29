join -1 2 -2 2 <(sort -k2,2 bulkbc.csv) <(sort -k2,2 Bulk_DNA_1_scars.txt) > Bulk_DNA_1_scars_bbc.txt
sort Bulk_DNA_1_scars_bbc.txt | uniq -c > Bulk_DNA_1_scar_bccounts.txt
rm Bulk_DNA_1_scars_bbc.txt
