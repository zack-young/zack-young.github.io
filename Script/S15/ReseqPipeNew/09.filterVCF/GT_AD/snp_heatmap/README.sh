#!/usr/bin/env bash

# one gene at a time
# -a is the metadata file. all runs need it.
# -g is the gene ID, and then follows group files
sh script.sh -a metadata_all.txt -g TraesCS1A01G000400 -n TraesCS1A01G000400 ../groupF/DF ../groupF/YC
sh script_region.sh -a metadata_all.append.txt -c chr6A -s 497000000 -e 563000000 use_samplelist &
# -f flanking region, including up and down stream
sh script.sh -a metadata_all.txt -g TraesCS1A01G000400 -f 2000 ../groupF/DF ../groupF/YC

# run 10 jobs in parallel
parallel -j 10 --resume-failed --retries 5 --joblog parallel.log sh script.sh -a metadata_all.txt -g {} ../groupF/DF ../groupF/YC :::: test_gene_list.txt > sub.e
