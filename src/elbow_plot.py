#######
# This scripts computes the VB lower bound for varying number of donors
#######

import vireoSNP
import numpy as np

from scipy import sparse
from scipy.io import mmread

# This avoids using X 
import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
from vireoSNP.plot.base_plot import heat_matrix

import argparse

# Read in Arguments
parser = argparse.ArgumentParser()
parser.add_argument('--refmatrix',  
                    help='Path to the ref matrix from vartrix')
parser.add_argument('--altmatrix', 
                    help='Path to the alt matrix from vartrix')
parser.add_argument('--cells', 
                    help='Path to barcodes containing cells')
parser.add_argument('--out', 
                    help='Path to png for plot')

args = parser.parse_args()
alt = args.altmatrix
ref = args.refmatrix
cells = args.cells
out_plot = args.out

# Load Vartrix data
vtrix = vireoSNP.read_vartrix(alt, ref, cells)
AD = vtrix['AD']
DP = vtrix['DP']

# Run vireo for varying number of donors
n_donor_list = np.arange(2, 7)
ELBO_list_all = []
for n_don in n_donor_list:
    res = vireoSNP.vireo_wrap(AD, DP, n_donor = n_don, learn_GT=True,
            n_extra_donor = 0, ASE_mode=False, fix_beta_sum=False,
            n_init=50, check_doublet=True, random_seed=1)
    ELBO_list_all.append(res['LB_list'])


# Plot
fig = plt.figure(figsize=(5, 4), dpi=100)
plt.plot(n_donor_list - 1, np.max(ELBO_list_all, axis=1))
plt.boxplot(ELBO_list_all)
plt.xticks(n_donor_list - 1, n_donor_list)
plt.ylabel("ELBO")
plt.xlabel("n_clones")
plt.savefig(out_plot)
