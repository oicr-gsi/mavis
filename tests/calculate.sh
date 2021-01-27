#!/bin/bash
cd $1

find .  -name "*\.zip" -exec unzip -q {} \; >/dev/null # unzip the results files

# - .tab files have deterministic names, non-deterministic contents (for some columns)
# - .json/.svg files have non-deterministic names, deterministic contents

# Therefore:
# - Check contents (not names) of .json/.svg files by comparing sorted lists of md5 checksums
# - Check invariant columns only for .tab files

echo ".tab files:"
find . -name "mavis_summary_all_WT.EPT0068.tab" | xargs cut -f 2,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54 | md5sum
find . -name "mavis_summary_WT.EPT0068_non-synonymous_coding_variants.tab" | xargs cut -f 2,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54 | md5sum
find . -name "filtered_pairs.tab" | xargs cut -f 2,4,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,56,57,58,59,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,108,109 | md5sum
echo ".json files:"
find . -name "*\.json" | xargs md5sum | cut -f 1 -d " " | sort
echo ".svg files:"
find . -name "*\.svg" | xargs md5sum | cut -f 1 -d " " | sort

