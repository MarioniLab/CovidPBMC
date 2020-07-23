#! /usr/bin/bash

## fish out ADT indicies from the Undetermined reads
python3 scripts/retrieve_adt.py  --index-fastq test_data/SLX-19160/fastq/Undetermined_S0_L002_I1_001.fastq.gz --read1-fastq  test_data/SLX-19160/fastq/Undetermined_S0_L002_R1_001.fastq.gz --read2-fastq test_data/SLX-19160/fastq/Undetermined_S0_L002_R2_001.fastq.gz  --index ATCACGAT --output-prefix test_data/SLX-19160/SIGAF7_ADT_S0_L002
