"""
Code to generate large amount of parton-level
information in MadGraph and pipe it through a
Pythia script.
"""

import sys, os
import random
import numpy as np

mg_dir = "/group/hepheno/smsharma/Dark-Showers/MG5_aMC_v2_5_6_patch/bin/" # MadGraph directory

n_jobs = 1 # Number of jobs to submit
prefix = "sig_tchannel_em" # Prefix that the MG job was generated with
nevt_each = 10000 # Number of events in each job
nmatch = 2 # Pythia nmatch parameter

batch='''#!/bin/bash
#SBATCH -N 1   # node count
#SBATCH --ntasks-per-node=12
#SBATCH --mem=48gb
#SBATCH -t 01:00:00
##SBATCH --mail-type=begin
##SBATCH --mail-type=end
##SBATCH --mail-user=smsharma@princeton.edu
#SBATCH -p hepheno
source /opt/rh/devtoolset-3/enable
source activate venv_py27
source /group/hepheno/smsharma/Dark-Showers/env.sh
source /group/hepheno/heptools/root/bin/thisroot.sh
cd '''

# for i in [10,20,50,100]:
for i in np.arange(600,3000,200):
# for i in [500,1000,2000,3000,4000,5000,7000,10000,20000,50000,100000]:
# for i in [500,1000,3000,5000,7000,10000,50000,100000]:
# for i in [7000,10000,20000,50000,100000]:
# for i in [5]:
	# Copy MG folder for each job
	os.system("cp -r " 
		+ mg_dir
		+ prefix + " "
		+ mg_dir
		+ prefix + "_" + str(i))
	batchn = batch + mg_dir + prefix + "_" + str(i) + "/bin/\n" # Go to appropriate MG folder
	seed = random.randrange(0,123120) # Random seed
	batchn += "sed -i 's/<phimass>/"+str(i)+"/g' ../Cards/param_card.dat\n" # Change seed in MG run card
	batchn += "sed -i 's/0   = iseed/"+str(seed)+"   = iseed/g' ../Cards/run_card.dat\n" # Change seed in MG run card
	batchn += "sed -i 's/10000 = nevents/"+str(nevt_each)+" = nevents/g' ../Cards/run_card.dat\n" # Change events in MG run card
	batchn += "rm -r " + mg_dir + prefix + "_" + str(i) + "/RunWeb\n" # Remove random dumb file
	batchn += "./generate_events -f\n" # Generate Pythia events
	batchn += "gunzip " + mg_dir + prefix+"_" + str(i) +"/Events/run_01/unweighted_events.lhe.gz\n" # Unzip created LHE file
	batchn += "sed -i 's/49001010/4900101/g' ../Events/run_01/unweighted_events.lhe\n" # Change events in MG run card	
	batchn += "sed -i 's/49001011/4900101/g' ../Events/run_01/unweighted_events.lhe\n" # Change events in MG run card	
	batchn += "sed -i 's/49001012/4900101/g' ../Events/run_01/unweighted_events.lhe\n" # Change events in MG run card	
	batchn += "sed -i 's/49001013/4900101/g' ../Events/run_01/unweighted_events.lhe\n" # Change events in MG run card	
	batchn += "sed -i 's/49001014/4900101/g' ../Events/run_01/unweighted_events.lhe\n" # Change events in MG run card	
	# batchn += "cd /group/hepheno/smsharma/Dark-Showers/gen/\n" # Go to Pythia script folder
	# batchn += "./monojet.exe -m lhe -w -i " + mg_dir + prefix + "_" + str(i) + "/Events/run_01/unweighted_events.lhe -metmin 0 -n 500000 -o monojet_mphi_" + str(i) + "-phimass " + str(i) + "-lambda " + str(i/2) + "\n"
	fname = "batch/batch_tch_" + prefix + "_" + str(i) + ".batch" # 
	f=open(fname, "w")
	f.write(batchn)
	f.close()
	os.system("sbatch " + fname);
