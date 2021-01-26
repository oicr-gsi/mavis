#!/bin/bash
cd $1

find .  -name "*\.zip" -exec unzip {} \; >/dev/null # unzip the results files

# filenames will vary, but contents should stay the same
# approach 1:
# find the md5 sums only and sort them

echo "Method 1:"
find . -type f -not -path "./*.zip" | xargs md5sum | cut -f 1 -d " " | sort

# approach 2:
# verify checksums of files with invariant names (eg. ./mavis_summary_all_WT.EPT0068.tab), only check for numbers of others (ie. total numbers of JSON/SVG files)

echo "Method 2:"
find . -name "\*.tab" | sort | xargs md5sum
find . -name "\*.json" | xargs wc -l
find . -name "\*.svg" | xargs wc -l

# run test overnight and pick one
