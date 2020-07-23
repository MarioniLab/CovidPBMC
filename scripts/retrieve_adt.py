#! /usr/bin/env python3

## This script finds read indexes matching the antibody library and generates the I1, R1 and R2 FASTQ files
import pysam
import os
import sys
import re
import logging
import argparse
import Bio.SeqIO as IO
import gzip

parser = argparse.ArgumentParser(description="Extract reads from ADT library given an index")

parser.add_argument("--index", dest="index_seq", type=str,
                    help="The sequence index used for the ADT library")

parser.add_argument("--index-fastq", dest="index_fastq", type=str,
                    help="Path to the FASTQ of indices")

parser.add_argument("--read1-fastq", dest="read1_fastq", type=str,
                    help="Path to the FASTQ for read1")

parser.add_argument("--read2-fastq", dest="read2_fastq", type=str,
                    help="Path to the FASTQ for read2")

parser.add_argument("--output-prefix", dest="out_prefix", type=str,
                    help="Prefix to use for output FASTQ files")

parser.add_argument("--log", dest="logfile", type=str,
                    help="Logging file destination", default=sys.stdout)

args = parser.parse_args()

if type(args.logfile) == str:
    logging.basicConfig(level=logging.INFO,
                        filename=args.logfile)
else:
    logging.basicConfig(level=logging.INFO)

index_regex = re.compile(args.index_seq)
match_reads = set()
index_reads = []
read1_reads = []
read2_reads = []

read_counter = 0
# reads are in the same order in each file, so I actually need the index of the reads
# I can also write the output FASTQ files on the fly
out_index = "{}_I1_001.fastq".format(args.out_prefix)
out_read1 = "{}_R1_001.fastq".format(args.out_prefix)
out_read2 = "{}_R2_001.fastq".format(args.out_prefix)

logging.info("Reading Index FASTQ file: {}".format(args.index_fastq))
with gzip.open(args.index_fastq, "rt") as fastq_handle:
    with open(out_index, "w") as ofile_index:
        for record in IO.parse(fastq_handle, "fastq"):
            title = record.description
            seq = str(record.seq)
            # find the read names that match the input index regex
            if index_regex.search(seq):
                match_reads.add(read_counter) # match_reads is a set for faster operations downstream
                # output matching reads on the fly to prevent excess memory consumption
                #logging.info("Writing FASTQ index record {} to file: {}".format(title, out_index))
                IO.write(record, ofile_index, "fastq")
                #index_reads.append(record)

            if read_counter % 1000000 == 0:
                logging.info("Read {} reads from Read1 file {}:".format(read_counter, args.index_fastq))
            read_counter += 1

logging.info("Found {} index reads matching in put index: {}".format(len(match_reads), args.index_seq))

#logging.info("Writing FASTQ index records to file: {}".format(out_index))
#print(index_reads)
#with open(out_index, "w") as ofile_index:
#    IO.write(index_reads, ofile_index, "fastq")

logging.info("Reading Read1 FASTQ file: {}".format(args.read1_fastq))
r1_counters = 0
# iterate over the other files and extract the relevant records
# storing all of these reads uses too much memory, so I'll output on the fly
with gzip.open(args.read1_fastq, "rt") as readq_handle:
    with open(out_read1, "w") as ofile_read1:
        for r1_record in IO.parse(readq_handle, "fastq"):
            r1_title = r1_record.description
            # set operations scale better for large lists
            if len(match_reads.intersection([r1_counters])):
                #logging.info("Writing FASTQ Read1 record {} to file: {}".format(r1_record.description, out_read1)) 
                IO.write(r1_record, ofile_read1, "fastq")
                #read1_reads.append(r1_record)

            if r1_counters % 1000000 == 0:
                logging.info("Read {} reads from Read1 file {}:".format(r1_counters, args.read1_fastq))
            r1_counters += 1

#logging.info("Writing {} FASTQ Read1 records to file: {}".format(len(read1_reads), out_read1))
#with open(out_read1, "w") as ofile_read1:
#    IO.write(read1_reads, ofile_read1, "fastq")

logging.info("Reading Read2 FASTQ file: {}".format(args.read2_fastq))
r2_counters = 0
with gzip.open(args.read2_fastq, "rt") as read2_handle:
    with open(out_read2, "w") as ofile_read2:
        for r2_record in IO.parse(read2_handle, "fastq"):
            r2_title = r2_record.description
            if len(match_reads.intersection([r2_counters])):
                #logging.info("Writing FASTQ Read2 record {} to file: {}".format(r2_record.description, out_read2))
                IO.write(r2_record, ofile_read2, "fastq")
                #read2_reads.append(r2_record)

            if r2_counters % 1000000 == 0:
                logging.info("Read {} reads from Read2 file {}:".format(r2_counters, args.read2_fastq))
            r2_counters += 1

#logging.info("Writing {} FASTQ Read2 records to file: {}".format(len(read2_reads), out_read2))
#with open(out_read2, "w") as ofile_read2:
#    IO.write(read2_reads, ofile_read2, "fastq")

