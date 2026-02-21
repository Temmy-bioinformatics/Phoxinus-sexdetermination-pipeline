# 08-extract-significant-kmers

# This is better ran on terminal and should be done ahead of running the final part of the snakefiles

# Find the value of the lowest p-value, the value printed is the most significant and should then be used as input in the next step 

cut -f 9 kmers.assoc.tab | grep -v 'P' | sort -g | head -1 # this takes some time

# Select these from the PLINK output
awk '$9 <= 3.414e-07' kmers.assoc.tab > most_significant_assoc.tab

# Most significant kmers are in the file 'most_significant_assoc.tab'.

# use this output table to create a presence/absence table of associated kmers for all samples

cat most_significant_assoc.tab | cut -f 2 > most_assoc_kmers.list
