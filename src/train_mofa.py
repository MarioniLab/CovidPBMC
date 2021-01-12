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
from scipy import io as scio

parser = argparse.ArgumentParser(description="Run MOFA on input multi-omics matrices")

parser.add_argument("--counts-list", dest="counts_files", type=str,
                    help="A comma-separated list of paths to input matrices (txt format)")

parser.add_argument("--omic", dest="omic", type=str,
                    choices=["ADT", "GEX", "joint"],
                    help="The type of input omic, ADT, GEX, or joint for both")

parser.add_argument("--donor-list", dest="donor_files", type=str,
                    help="A comma-separated list of paths to donor ID information, "
                    "must be of the same length as the number of files")

parser.add_argument("--gex-rownames", dest="gex_rownames", type=str,
                    default=None, help="A path to a file containing the GEX rownames")

parser.add_argument("--adt-rownames", dest="adt_rownames", type=str,
                    default=None, help="A path to a file containing the ADT rownames")

parser.add_argument("--gex-colnames", dest="gex_colnames", type=str,
                    default=None, help="A path to a file containing the GEX colnames")

parser.add_argument("--adt-colnames", dest="adt_colnames", type=str,
                    default=None, help="A path to a file containing the ADT colnames")

parser.add_argument("--format", dest="data_format", type=str,
                    choices=["wide", "long"], default="wide",
                    help="Input data format, either multiple data matrices in wide format, or a "
                    "single long format data frame.")

parser.add_argument("--group-column", dest="group_col", type=str,
                    help="Column from donors file that will be used for MOFA groups")

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

logging.info("Setting up numpy multi-threading. Using {} threads".format(args.use_threads))
os.environ['OPENBLAS_NUM_THREADS'] = str(args.use_threads)

logging.info("Reading in ADT counts matrices")

# initial the model entry point
mofa_ent = entry_point()

# each matrix should be gzipped
file_list = args.counts_files.split(",")

# read in the donor data frames as dict for reference
# must be the same number of files unless running joint mode
logging.info("Readining in donor information: {}".format(args.donor_files))
donor_list = args.donor_files.split(",")

if len(donor_list) != len(file_list):
    if args.omic == 'joint' and len(donor_list) == 1:
        donor_names = ["Gene", "ADT"]
        donor_df_list = [pd.read_table(D, sep="\t", header=0, index_col='CellID') for D in donor_list]
        donor_dict = dict(zip(donor_names,
                              donor_df_list * len(donor_names)))
    else:
        raise("The number of donor files is not the same as the input files")
else:
    donor_names = [FX.split("/")[-2] for FX in donor_list]
    donor_dict = dict(zip(donor_names,
                          [pd.read_table(D, sep="\t", header=0, index_col=None) for D in donor_list]))

# keep only the singlets
keep_dict = {}

# the best way to run this is to melt the whole data matrix into long format and combine with the 
# relevant meta data

matrix_list = []

# read in dimension names as pandas series
# need to change ".1" to "-1" artefact from R
logging.info("Reading in data dimension names")
dim_dict = {}
logging.info("GEX features: {}".format(args.gex_rownames))
gex_rows = pd.read_table(args.gex_rownames, sep="\t", index_col=None, header=None)
logging.info("GEX cell IDs: {}".format(args.gex_colnames))
gex_cols = pd.read_table(args.gex_colnames, sep="\t", index_col=None, header=None)
logging.info("ADT features: {}".format(args.adt_rownames))
adt_rows = pd.read_table(args.adt_rownames, sep="\t", index_col=None, header=None)
logging.info("ADT cell IDs: {}".format(args.adt_colnames))
adt_cols = pd.read_table(args.adt_colnames, sep="\t", index_col=None, header=None)

dim_dict['ADT'] = {'rows': adt_rows.iloc[:, 0].values,
                   'cols': [ax.replace('.1', '-1') for ax in adt_cols.iloc[:, 0].values]}
dim_dict['GEX'] = {'rows': gex_rows.iloc[:, 0].values,
                   'cols': [gx.replace('.1', '-1') for gx in gex_cols.iloc[:, 0].values]}

logging.info("Populating dimension dictionary, keys: {}".format(dim_dict.keys()))

for x in range(len(file_list)):
    x_file = file_list[x]
    # extract the sample name from the matrix file
    # assumes the first column contains the relevant feature IDs
    if args.data_format == "wide":
        if re.search("gz$", x_file):
            x_matrix = pd.read_table(x_file, delimiter="\t", compression="gzip",
                                     header=0, index_col=None)
        elif re.search("mtx$", x_file):
            logging.info("Reading matrix file: {}".format(x_file))
            x_matrix = scio.mmread(x_file)
            if args.omic == "joint":
                if re.search("ADT", x_file):
                    logging.info("Adding ADT dimension names: {}".format(x_matrix.shape))
                    x_matrix = pd.DataFrame.sparse.from_spmatrix(x_matrix,
                                                                 index=dim_dict['ADT']['rows'],
                                                                 columns=dim_dict['ADT']['cols'])
                    x_matrix["ADT"] = dim_dict['ADT']['rows']

                elif re.search("GEX", x_file):
                    logging.info("Adding GEX dimension names: {}".format(x_matrix.shape))
                    x_matrix = pd.DataFrame.sparse.from_spmatrix(x_matrix,
                                                                 index=dim_dict['GEX']['rows'],
                                                                 columns=dim_dict['GEX']['cols'])
                    x_matrix["Gene"] = dim_dict['GEX']['rows']

            elif args.omic == "GEX":
                logging.info("Adding GEX dimension names: {}".format(x_matrix.shape))
                x_matrix = pd.DataFrame.sparse.from_spmatrix(x_matrix,
                                                             index=dim_dict['GEX']['rows'],
                                                             columns=dim_dict['GEX']['cols'])
                x_matrix["Gene"] = dim_dict['GEX']['rows']

            elif args.omic == "ADT":
                logging.info("Adding ADT dimension names: {}".format(x_matrix.shape))
                x_matrix = pd.DataFrame.sparse.from_spmatrix(x_matrix,
                                                             index=dim_dict['ADT']['rows'],
                                                             columns=dim_dict['ADT']['cols'])
                x_matrix["ADT"] = dim_dict['ADT']['rows']
        else:
            x_matrix = pd.read_table(x_file, delimiter="\t",
                                     header=0, index_col=None)

        logging.info("Melting matrices into long data frames")
        # melt based on omic type
        if args.omic == "ADT":
            melt_id = "ADT"
        elif args.omic == "GEX":
            melt_id = "Gene"
        elif args.omic == "joint":
            # infer from filename
            if re.search("ADT", x_file):
                melt_id = "ADT"
            elif re.search("GEX", x_file):
                melt_id = "Gene"

        # need to test if key is in the columns, if not then use the index
        try:
            x_melt = x_matrix.melt(id_vars=melt_id)
        except KeyError:
            x_melt = x_matrix.melt()
        print(x_melt.head())

        x_melt.columns = ["feature", "sample", "value"]
        x_melt.loc[:, "view"] = melt_id
        if args.omic == "joint":
            donor_df = donor_dict[melt_id]
            # remove cell Ids from the donor list that don't appear in the ADT data
            donor_df = donor_df.loc[donor_df.index.isin(x_melt["sample"].unique()), :]
            x_melt = x_melt.merge(donor_df.loc[:, [args.group_col]], how='left', left_on='sample', right_index=True)
            x_melt.columns = ["feature", "sample", "value", "view", "group"]
        else:
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

    # check there are no NaNs in the groups - this has caused a massive headache.
    print(dt_frame.loc[dt_frame["group"].isnull(), "sample"].str.extract(pat="(BGCV[0-9]+)_\S+"))
    dt_frame.loc[dt_frame["group"].isnull(), "group"] = dt_frame.loc[dt_frame["group"].isnull(), "sample"].str.extract(pat="(BGCV[0-9]+)_\S+").values
    print(dt_frame[dt_frame["group"].isnull()].head())

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
print(dt_frame.head())
mofa_ent.set_data_df(dt_frame)
#mofa_ent.set_data_df(dt_frame, likelihoods = ["gaussian" for x in range(len(file_list))])

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
    startELBO = 1, 
    freqELBO = 1, 
    dropR2 = None, 
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
