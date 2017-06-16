"""
Code to generate large amount of parton-level
information in MadGraph and pipe it through a
Pythia script.
"""

import sys, os
import random

base_dir = "/group/hepheno/smsharma/Dark-Showers/"

mg_dir = base_dir + "/MG5_aMC_v2_5_6_patch/bin/" # MadGraph directory
data_dir = base_dir + "/data/"

# Which backgrounds should generate (gridpacks should exist in these dirs already)
# bkg_prefix = ["bkg_Zj", "bkg_Wj", "bkg_ttbar_semilep", "bkg_ttbar_dilep"]
bkg_prefix = ["bkg_qcd", "bkg_ttbar_dilep"]


batch='''#!/bin/bash
#SBATCH -N 1   # node count
#SBATCH --ntasks-per-node=12
#SBATCH --mem=36gb
#SBATCH -t 01:00:00
##SBATCH --mail-type=begin
##SBATCH --mail-type=end
##SBATCH --mail-user=smsharma@princeton.edu
#SBATCH -p dept
source /opt/rh/devtoolset-3/enable
source activate venv_py27
source /group/hepheno/smsharma/Dark-Showers/env.sh
source /group/hepheno/heptools/root/bin/thisroot.sh
cd '''

for ibkg, prefix in enumerate(bkg_prefix): 
			
	batchn = batch + mg_dir + prefix + "\n" # Go to appropriate MG folder
	batchn += "rm RunWeb" +"\n" # Delete random shit file if needed			

	# Copy over gridpack to temp
	batchn += "./bin/generate_events -f" +"\n" # Clean up previous gridpack
	batchn += "tar -xvf run_01_gridpack.tar.gz" + "\n" # Clean up previous gridpack
	batchn += "./madevent/bin/compile" +"\n" # Make folder for this run
	
	fname = "batch/batch_gridpack_" + prefix + "_.batch" # 
	f=open(fname, "w")
	f.write(batchn)
	f.close()
	os.system("sbatch " + fname);