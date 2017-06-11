"""
Code to generate large amount of parton-level
information in MadGraph and pipe it through a
Pythia script.
"""

import sys, os
import random

import numpy as np

mg_dir = "/group/hepheno/smsharma/Dark-Showers/MG5_aMC_v2_5_2/bin/" # MadGraph directory

n_jobs = 1 # Number of jobs to submit
prefix = "DM_Zprime" # Prefix that the MG job was generated with
nevt_each = 50000 # Number of events in each job
nmatch = 2 # Pythia nmatch parameter

batch='''#!/bin/bash
#SBATCH -N 1   # node count
#SBATCH --ntasks-per-node=16
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

# for MZ in [500,1000,3000]:
# for MZ in [5000, 8000, 9000]:
# for MZ in [200,500,1000]:

gX = 1
gq = 0.1
# masses = np.arange(1000,4000,100)
masses = [4000]
for MZ in masses:
	# Copy MG folder for each job
	os.system("cp -r " 
		+ mg_dir
		+ prefix + " "
		+ mg_dir
		+ prefix + "_" + str(MZ))
	batchn = batch + mg_dir + prefix + "_" + str(MZ) + "/bin/\n" # Go to appropriate MG folder
	seed = random.randrange(0,123120) # Random seed
	print "The seed is", seed
	batchn += "sed -i 's/1.000000e+03 # MY1/"+str(MZ)+" # MY1/g' ../Cards/param_card.dat\n" # Change seed in MG run card
	W = MZ/3*(2*2*gX**2 + 6*3*gq**2)/(4*np.pi) # Calculate width
	batchn += "sed -i 's/1.000000e+01 # WY1/"+str(W)+" # WY1/g' ../Cards/param_card.dat\n" # Change seed in MG run card
	batchn += "sed -i 's/0 = iseed/"+str(seed)+" = iseed/g' ../Cards/run_card.dat\n" # Change seed in MG run card
	batchn += "sed -i 's/100000 = nevents/"+str(nevt_each)+" = nevents/g' ../Cards/run_card.dat\n" # Change events in MG run card
	batchn += "rm -r " + mg_dir + prefix + "_" + str(MZ) + "/RunWeb\n" # Remove random dumb file
	batchn += "./generate_events -f\n" # Generate Pythia events
	batchn += "gunzip " + mg_dir + prefix+"_" + str(MZ) +"/Events/run_01/unweighted_events.lhe.gz\n" # Unzip created LHE file
	batchn += "sed -i 's/5000521/4900101/g' ../Events/run_01/unweighted_events.lhe\n" # Change events in MG run card	
	batchn += "cd /group/hepheno/smsharma/Dark-Showers/gen/\n" # Go to Pythia script folder
	# batchn += "./monojet.exe -m lhe -w -MZ " + mg_dir + prefix + "_" + str(MZ) + "/Events/run_01/unweighted_events.lhe -metmin 0 -n 500000 -o monojet_mphi_" + str(MZ) + "-phimass " + str(MZ) + "-lambda " + str(MZ/2) + "\n"
	fname = "batch/batch_" + prefix + "_" + str(MZ) + ".batch" # 
	f=open(fname, "w")
	f.write(batchn)
	f.close()
	os.system("sbatch " + fname);
