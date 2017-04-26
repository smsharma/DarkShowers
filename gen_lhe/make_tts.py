"""
Code to generate large amount of parton-level
information in MadGraph and pipe it through a
Pythia script.
"""

import sys, os
import random

mg_dir = "/group/hepheno/smsharma/Dark-Showers/MG5_aMC_v2_5_2/bin/" # MadGraph directory

n_jobs = 1 # Number of jobs to submit
prefix = "bkg_ttbar_semilep" # Prefix that the MG job was generated with
nevt_each = 100000 # Number of events in each job
nmatch = 5 # Pythia nmatch parameter

batch='''#!/bin/bash
#SBATCH -N 1   # node count
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4gb
#SBATCH -t 16:00:00
##SBATCH --mail-type=begin
##SBATCH --mail-type=end
##SBATCH --mail-user=smsharma@princeton.edu
#SBATCH -p hepheno
source /opt/rh/devtoolset-3/enable
source activate venv_py27
source /group/hepheno/smsharma/Dark-Showers/env.sh
source /group/hepheno/heptools/root/bin/thisroot.sh
cd '''

for i in range(n_jobs):
	# Copy MG folder for each job
	os.system("cp -r " 
		+ mg_dir
		+ prefix + " "
		+ mg_dir
		+ prefix + "_" + str(i))
	batchn = batch + mg_dir + prefix + "_" + str(i) + "/bin/\n" # Go to appropriate MG folder
	seed = random.randrange(0,123120) # Random seed
	batchn += "sed -i 's/0   = iseed/"+str(seed)+"   = iseed/g' ../Cards/run_card.dat\n" # Change seed in MG run card
	batchn += "sed -i 's/50000 = nevents/"+str(nevt_each)+" = nevents/g' ../Cards/run_card.dat\n" # Change events in MG run card
	batchn += "rm -r " + mg_dir + prefix + "_" + str(i) + "/RunWeb\n" # Remove random dumb file
	batchn += "./generate_events -f\n" # Generate Pythia events
	batchn += "gunzip " + mg_dir + prefix+"_" + str(i) +"/Events/run_01/unweighted_events.lhe.gz\n" # Unzip created LHE file
	batchn += "cd /group/hepheno/smsharma/Dark-Showers/gen/\n" # Go to Pythia script folder
	batchn += "./monojet.exe -m lhe -w -i " + mg_dir + prefix + "_" + str(i) + "/Events/run_01/unweighted_events.lhe -metmin 0 -Zprime -n 500000 -o bkg_ttbar_semilep_" + str(i) + " -nmatch " + str(nmatch) + "\n"
	fname = "batch/batch_" + prefix + "_" + str(i) + ".batch" # 
	f=open(fname, "w")
	f.write(batchn)
	f.close()
	os.system("sbatch " + fname);