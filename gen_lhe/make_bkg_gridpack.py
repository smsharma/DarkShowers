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
# bkg_prefix = ["bkg_qcd", "bkg_Zj", "bkg_Wj", "bkg_ttbar_semilep", "bkg_ttbar_dilep"]
bkg_prefix = ["bkg_qcd"]

nevents = 5000 # Events per job
ijobs = [200, 0, 0, 0] # Start index
njobs = [400, 0, 0,0] # Number of jobs to submit

# nmatch = [4, 2, 2, 5, 3] # Pythia nmatch parameter
nmatch = [4] # Pythia nmatch parameter
recluster = False # Whether to recluster jets in pythia

if recluster:
	recluster_str = " -Zprime "
else:
	recluster_str = " "

batch='''#!/bin/bash
#SBATCH -N 1   # node count
#SBATCH --ntasks-per-node=1
#SBATCH --mem=3gb
#SBATCH -t 02:00:00
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

	print "Submitting to generate", nevents*njobs[ibkg], "events for", prefix, "! Nice"

	# Make folder for final ntuples
	if not os.path.exists(data_dir + prefix):
		os.mkdir(data_dir + prefix)

	seeds = random.sample(range(1, 30081), njobs[ibkg])

	for ij in range(ijobs[ibkg], ijobs[ibkg] + njobs[ibkg]):

		batchn = batch + mg_dir + prefix + "\n" # Go to appropriate MG folder

		# Copy over gridpack to temp
		batchn += "rm -rf /tmp/sid_tmp/gridpack_" +  str(ij) +"\n" # Clean up previous gridpack
		batchn += "mkdir /tmp/sid_tmp/" +"\n" # Make folder for this run
		batchn += "mkdir /tmp/sid_tmp/gridpack_" +  str(ij) +"\n" # Make folder for this run
		batchn += "cp -r run.sh madevent /tmp/sid_tmp/gridpack_"  + str(ij) +"\n" # Copy over gridpack files
		batchn += "cd /tmp/sid_tmp/gridpack_" + str(ij) +"\n" # Go to folder
		batchn += "rm madevent/RunWeb" +"\n" # Delete random shit file if needed			
		
		# Generate MG events
		batchn += "./run.sh " + str(nevents) + " " + str(seeds[ij - ijobs[ibkg]]) +"\n" 
		
		# Munge LHE file
		lhe_file = "/tmp/sid_tmp/gridpack_" + str(ij) + "/events.lhe"
		batchn += "gunzip " + lhe_file + ".gz\n" #
		batchn += "cp " + lhe_file + " " + mg_dir + prefix + "/Events/events_" + str(ij) + ".lhe" "\n" #
		lhe_file = mg_dir + prefix + "/Events/events_"+ str(ij) + ".lhe"

		# Clean up
		batchn += "rm -rf /tmp/sid_tmp/gridpack_" +  str(ij) + "\n" 		

		# Run pythia script
		batchn += "cd " + base_dir + "/gen/\n" # Go to Pythia script folder
		batchn += "./monojet.exe -m lhe -w -i " + lhe_file + " -metmin 0 " + recluster_str + " -n " + str(nevents) + " -o " + data_dir + prefix + "/events_" + str(ij) + " -nmatch " + str(nmatch[ibkg]) + "\n"
		
		fname = "batch/batch_" + prefix + "_" + str(ij) + ".batch" # 
		f=open(fname, "w")
		f.write(batchn)
		f.close()
		os.system("sbatch " + fname);