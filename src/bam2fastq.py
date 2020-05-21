#! /usr/bin/env Python3
# convert input bam files to FASTQ using pysam

import os
import sys
import pysam
import argparse
import logging
import gzip


def parse_bam(bam, bam_filter=False):
    '''
    Parse the input bam file and return a generator of sequence alignments for 
    just in time evaluation.
    '''

    align = pysam.AlignmentFile(bam, "rb")
    for record in align.fetch(until_eof=True):

        yield(record)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Convert a BAM file to FASTQ format")

    parser.add_argument("--bam", dest="bam_file", type=str,
                        help="Input bam file")

    parser.add_argument("--fastq", dest="fastq_file", type=str,
                        help="Output FASTQ file")

    parser.add_argument("--filter-pairs", dest="filter_pair", action='store_true',
                        help="Flag to keep paired reads only")

    parser.add_argument("--log", dest="logfile", type=str, default=sys.stdout,
                        help="Path to logging file")

    args = parser.parse_args()

    # setup the logger
    if type(args.logfile) == str:
        logging.basicConfig(level=logging.INFO,
                            filename=args.logfile)
    else:
        logging.basicConfig(level=logging.INFO)
                        
    logging.info("Parsing bam file: {}".format(args.bam_file))
    parse_records = parse_bam(args.bam_file)

    reads_counted = 0
    with gzip.open(args.fastq_file, "wb") as fqfile:
        logging.info("Only retaining proper reads")
        for x in parse_records:
            # this conversion will need to keep the cell barcode, UMI and
            # index in the final read name for later
            # format <readname>.<cellbarcode>.<umi>.<index>
            # CB and UB are the error-corrected cell barcode and UMI, respectively
            # check for read tags
            if x.has_tag('RG'):
                readname = x.get_tag('RG')
            else:
                readname = False
        
            if x.has_tag('CB'):
                barcode = x.get_tag('CB')
            else:
                barcode = False

            if x.has_tag('UB'):
                umi = x.get_tag('UB')
            else:
                umi = False

            if x.has_tag('BC'):
                index = x.get_tag('BC')
            else:
                index = False

            if all([readname, barcode, umi, index]):
                # print out to FASTQ file on the fly
                fastq_name = "{}.{}.{}.{}".format(readname, barcode, umi, index)
                out_string = "{}\n{}\n+\n{}\n".format(fastq_name, x.seq, ''.join(map(lambda Q: chr( Q+33), x.query_qualities)))
                fqfile.write(out_string.encode())
            else:
                if (reads_counted % 10000) == 0:
                    reads_counted += 1
                    logging.warning("{} reads rejected for missing information".format(reads_counted))
            
