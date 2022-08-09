#!/bin/bash
cd $1

# - Check contents (not names) of .json/.svg files by comparing sorted lists of md5 checksums
# - Check invariant columns only for .tab files

echo ".tab files:"
find . -name "*.mavis_summary.tab" | xargs cut -f 2,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54 | md5sum
find . -name "*.WT_non-synonymous_coding_variants.tab" | xargs cut -f 2,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54 | md5sum


find .  -name "*\.zip" -exec unzip -q {} \; >/dev/null # unzip the results files
echo ".json files:"
find . -name "*\.json" | xargs md5sum | cut -f 1 -d " " | sort
echo ".svg files:"
find . -name "*\.svg" | xargs md5sum | cut -f 1 -d " " | sort
