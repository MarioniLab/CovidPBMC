#! /usr/bin/env python3

## Run MOFA on either a single sample or multiple samples, each of which is treated as a separate group

import os
import sys
import re
import logging
import argparse
import gzip
import mofapy2
from mofapy2.run.entry_point import entry_point
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(description="Extract reads from ADT library given an index")

parser.add_argument("--counts-list", dest="counts_files", type=str,
                    help="A comma-separated list of paths to ADT counts matrices (txt format)")

parser.add_argument("--donor-list", dest="donor_files", type=str,
                    help="A comma-separated list of paths to donor ID information")

parser.add_argument("--format", dest="data_format", type=str,
                    choices=["wide", "long"], default="wide",
                    help="Input data format, either multiple data matrices in wide format, or a "
                    "single long format data frame.")

parser.add_argument("--factors", dest="factors", type=int,
                    help="The number of factors to learn")

parser.add_argument("--ignore-groups", dest="ignore_groups", action="store_true",
                    default=False, help="If present, ignore the group column and run all "
                    "samples as a single group")

parser.add_argument("--runmode", dest="runmode", type=str,
                    choices=["fast", "medium", "slow"],
                    default="medium",
                    help="Mode to run MOFA in. Either fast, medium or slow")

parser.add_argument("--MOFAoutput", dest="MOFA_out", type=str,
                    help="Path to HDF5 file to save trained MOFA model")

parser.add_argument("--output-prefix", dest="out_prefix", type=str,
                    help="Prefix to use for output files of denoised expression matrices")

parser.add_argument("--log", dest="logfile", type=str,
                    help="Logging file destination", default=sys.stdout)

parser.add_argument("--threads", dest="use_threads", type=int,
                    default=1,
                    help="The number of threads to use for parallelised MOFA")

args = parser.parse_args()

if type(args.logfile) == str:
    logging.basicConfig(level=logging.INFO,
                        filename=args.logfile)
else:
    logging.basicConfig(level=logging.INFO)

logging("Setting up numpy multi-threading. Using {} threads".format(args.use_threads))
os.environ['OPENBLAS_NUM_THREADS'] = str(args.use_threads)

logging.info("Reading in ADT counts matrices")

# initial the model entry point
mofa_ent = entry_point()

# each matrix should be gzipped
file_list = args.counts_files.split(",")

# read in the donor data frames as dict for reference
donor_list = args.donor_files.split(",")
donor_names = [FX.split("/")[-2] for FX in donor_list]

donor_dict = dict(zip(donor_names,
                      [pd.read_table(D, sep="\t", header=0, index_col=None) for D in donor_list]))

# keep only the singlets
keep_dict = {}

# the best way to run this is to melt the whole data matrix into long format and combine with the 
# relevant meta data

matrix_list = []

for x in range(len(file_list)):
    x_file = file_list[x]
    # extract the sample name from the matrix file
    # assumes the first column contains the relevant feature IDs
    if args.data_format == "wide":
        x_matrix = pd.read_table(x_file, delimiter="\t", compression="gzip",
                                 header=0, index_col=None)

        logging.info("Melting matrices into long data frames")
        x_melt = x_matrix.melt(id_vars="ADT")
        x_melt.columns = ["feature", "sample", "value"]
        x_melt.loc[:, "view"] = "ADT"
        x_melt.loc[:, "group"] = donor_names[x]
        print(x_melt.head())

    elif args.data_format == "long":
         x_melt = pd.read_table(x_file, delimiter="\t", compression="gzip",
                                 header=0, index_col=None)
         x_melt.columns = ["feature", "sample", "value", "group", "view"]
    else:
        raise AttributeError("Data format not recognised, must be either wide or long")

    matrix_list.append(x_melt)

    # make one big data frame
    dt_frame = pd.concat(matrix_list)

if args.ignore_groups:
    logging.info("Ignore groups label from data - running in single group mode")
    # set the group value to all 0
    dt_frame.group = "Group1"

# setup data options
mofa_ent.set_data_options(
    scale_groups = True, 
    scale_views = True
)

# add the data to the model
logging.info("Setting up MOFA data options")
mofa_ent.set_data_df(dt_frame, likelihoods = ["gaussian"])

logging.info("Setting up MOFA model options")
mofa_ent.set_model_options(
    factors = args.factors, 
    spikeslab_weights = True, 
    ard_factors = True,
    ard_weights = True
)

logging.info("Setting up training options")
mofa_ent.set_train_options(
    iter = 1000, 
    convergence_mode = args.runmode, 
    #drop_factor_threshold = -1,
    startELBO = 1, 
    freqELBO = 1, 
    dropR2 = 0.0000001, 
    gpu_mode = False,
    verbose = False,
    seed = 42
)

logging.info("Building MOFA model")
mofa_ent.build()

logging.info("Training MOFA model")
mofa_ent.run()

logging.info("Saving MOFA model to: {}".format(args.MOFA_out))
mofa_ent.save(args.MOFA_out)
