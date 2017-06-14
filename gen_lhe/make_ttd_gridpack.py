"""
Code to generate large amount of parton-level
information in MadGraph and pipe it through a
Pythia script.
"""

import sys, os
import random

base_dir = "/group/hepheno/smsharma/Dark-Showers/"

mg_dir = base_dir + "/MG5_aMC_v2_5_5/bin/" # MadGraph directory
data_dir = base_dir + "/data/"

nevents = 5000
njobs = 200 # Number of jobs to submit
prefix = "bkg_ttbar_dilep" # Prefix that the MG job was generated with
nmatch = 3 # Pythia nmatch parameter
recluster = False # Whether to recluster jets in pythia

if recluster:
	recluster_str = " -Zprime "
else:
	recluster_str = " "

print "Submitting to generate", nevents*njobs, "events! Nice"

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

os.system("rm -r " + mg_dir + prefix + "/gridpack_*") # Remove all previous runs

# Make folder for final ntuples
if not os.path.exists(data_dir + prefix):
	os.mkdir(data_dir + prefix)

for i in range(njobs):
	
	batchn = batch + mg_dir + prefix + "\n" # Go to appropriate MG folder
	batchn += "mkdir gridpack_" + str(seeds[i]) +"\n" # Make folder for this run
	batchn += "cp -r ./run.sh madevent gridpack_" + str(seeds[i]) +"\n" # Copy over gridpack files
	batchn += "cd gridpack_" + str(seeds[i]) +"\n" # Go to folder
	batchn += "rm madevent/RunWeb" +"\n" # Delete random shit file
	
	batchn += "./run.sh " + str(nevents) + " " + str(seeds[i]) +"\n" # Generate Pythia events
	
	lhe_file = mg_dir + prefix + "/gridpack_" + str(seeds[i]) + "/events.lhe"
	batchn += "gunzip " + lhe_file + ".gz\n" #
	
	batchn += "cd " + base_dir + "/gen/\n" # Go to Pythia script folder

	# Run pythia script
	batchn += "./monojet.exe -m lhe -w -i " + lhe_file + " -metmin 0 " + recluster_str + " -n " + str(nevents) + " -o " + data_dir + prefix + "/events_" + str(i) + " -nmatch " + str(nmatch) + "\n"
	
	fname = "batch/batch_" + prefix + "_" + str(i) + ".batch" # 
	f=open(fname, "w")
	f.write(batchn)
	f.close()
	os.system("sbatch " + fname);