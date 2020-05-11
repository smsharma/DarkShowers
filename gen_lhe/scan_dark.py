"""
Code to generate large amount of parton-level
information in MadGraph and pipe it through a
Pythia script.
"""

import sys, os
import random
import numpy as np

base_dir = "/group/hepheno/smsharma/Dark-Showers/"

mg_dir = base_dir + "/MG5_aMC_v2_5_6_patch/bin/" # MadGraph directory
data_dir = base_dir + "/data/"

# Which backgrounds should generate (gridpacks should exist in these dirs already)
# bkg_prefix = ["bkg_Zj", "bkg_Wj", "bkg_ttbar_semilep", "bkg_ttbar_dilep"]
bkg_prefix = ["bkg_qcd", "bkg_ttbar_dilep"]


batch='''#!/bin/bash
#SBATCH -N 1   # node count
#SBATCH --ntasks-per-node=1
#SBATCH --mem=3gb
#SBATCH -t 00:10:00
##SBATCH --mail-type=begin
##SBATCH --mail-type=end
##SBATCH --mail-user=smsharma@princeton.edu
#SBATCH -p dept
source /opt/rh/devtoolset-3/enable
source activate venv_py27
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/group/hepheno/heptools/HepMC-2.06.09/lib/
source /group/hepheno/smsharma/Dark-Showers/env.sh
source /group/hepheno/heptools/root/bin/thisroot.sh
cd '''

# lambda_range = np.linspace(0.01, 10, 20)
lambda_range = np.linspace(1, 400, 10)
rinv_range = np.linspace(0,1, 10)

for lambdo in lambda_range:
	for rinv in rinv_range:

			
		batchn = batch + base_dir +"/gen" + "\n" # Go to appropriate MG folder

		# Copy over gridpack to temp
		batchn += "./monojet.exe -m lhe -w -i /group/hepheno/smsharma/Dark-Showers//MG5_aMC_v2_5_6_patch_signal/bin/sig_zprime/Events/events_100_100000_0.lhe -o out_100_" + str(lambdo) + "_"+ str(rinv)[:5] +  " -n 1000 -metmin 0 -nmatch 2 -lambda " +str(lambdo) + " -inv " + str(rinv) +"\n" # Clean up previous gridpack
		
		fname = "batch/batch_gridpack_" + str(lambdo) + "_"+ str(rinv) + ".batch" # 
		f=open(fname, "w")
		f.write(batchn)
		f.close()
		os.system("sbatch " + fname);