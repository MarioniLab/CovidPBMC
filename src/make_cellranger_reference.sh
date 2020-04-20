#! /usr/bin/bash

## Make the cellranger format reference genome and annotations
## cellranger needs to be in the PATH variable
## 1) Genome ID
## 2) Path to Genome FASTA (uncompressed)
## 3) Path to Genome annotation GTF (uncompressed)
## 4) Optional reference string

cellranger mkref --genome=$1 \
                 --fasta=$2 \
                 --genes=$3 \
                 --ref-version=$4
