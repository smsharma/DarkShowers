"""
Code to generate large amount of parton-level
information in MadGraph and pipe it through a
Pythia script.
"""

import sys, os
import random

mg_dir = "/group/hepheno/smsharma/Dark-Showers/MG5_aMC_v2_5_2/bin/" # MadGraph directory

nevents = 5000
njobs = 30 # Number of jobs to submit
prefix = "bkg_Zj" # Prefix that the MG job was generated with
nmatch = 2 # Pythia nmatch parameter

batch='''#!/bin/bash
#SBATCH -N 1   # node count
#SBATCH --ntasks-per-node=1
#SBATCH --mem=3gb
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

seeds = random.sample(range(njobs), njobs)

for i in range(njobs):
	batchn = batch + mg_dir + prefix + "\n" # Go to appropriate MG folder
	batchn += "mkdir gridpack_" + str(seeds[i]) +"\n" # Generate Pythia events
	batchn += "cp -r ./run.sh madevent gridpack_" + str(seeds[i]) +"\n" # Generate Pythia events
	batchn += "cd gridpack_" + str(seeds[i]) +"\n" # Generate Pythia events
	batchn += "rm madevent/RunWeb" +"\n" # Generate Pythia events
	batchn += "./run.sh " + str(nevents) + " " + str(seeds[i]) +"\n" # Generate Pythia events
	# batchn += "gunzip " + mg_dir + prefix+"_" + str(i) +"/Events/run_01/unweighted_events.lhe.gz\n" # Unzip created LHE file
	# batchn += "cd /group/hepheno/smsharma/Dark-Showers/gen/\n" # Go to Pythia script folder
	# batchn += "./monojet.exe -m lhe -w -i " + mg_dir + prefix + "_" + str(i) + "/Events/run_01/unweighted_events.lhe -metmin 0 -Zprime -n 100000 -o bkg_Zj_" + str(i) + " -nmatch " + str(nmatch) + "\n"
	fname = "batch/batch_" + prefix + "_" + str(i) + ".batch" # 
	f=open(fname, "w")
	f.write(batchn)
	f.close()
	os.system("sbatch " + fname);