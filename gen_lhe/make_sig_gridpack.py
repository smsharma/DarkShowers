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

# Which signals should generate (gridpacks should exist in these dirs already)
bkg_prefix = ["sig_tchannel"]

nevents = 5000 # Events per job
ijobs = [0, 0, 0, 0, 0] # Start index
njobs = [4, 0, 0, 0, 0] # Number of jobs to submit
 # Prefix that the MG job was generated with
nmatch = [2, 2, 2, 5, 3] # Pythia nmatch parameter
recluster = False # Whether to recluster jets in pythia

masses = np.arange(500,3100,200)
rinvs = [0.0, 0.01,0.1,.2,.3,.4,.5,.6,.7,.8,.9, 0.99, 1.0]

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

	for imass, mass in enumerate(masses):

		seeds = random.sample(range(1, 30081), njobs[ibkg])

		for ij in range(ijobs[ibkg], ijobs[ibkg] + njobs[ibkg]):
			
			batchn = batch + mg_dir + prefix + "\n" # Go to appropriate MG folder
			
			# Copy over gridpack to temp
			batchn += "rm -rf /tmp/sid_tmp/gridpack_" + str(mass) + "_" +  str(ij) +"\n" # Clean up previous gridpack
			batchn += "mkdir /tmp/sid_tmp/" +"\n" # Make folder for this run
			batchn += "mkdir /tmp/sid_tmp/gridpack_"+ str(mass) + "_" +  str(ij) +"\n" # Make folder for this run
			batchn += "cp -r run.sh madevent /tmp/sid_tmp/gridpack_"+ str(mass) + "_"  + str(ij) +"\n" # Copy over gridpack files
			batchn += "cd /tmp/sid_tmp/gridpack_" + str(mass) + "_" + str(ij) +"\n" # Go to folder
			batchn += "rm madevent/RunWeb" +"\n" # Delete random shit file if needed			
			batchn += "sed -i 's/<PhiMass>/"+str(mass)+"/g' madevent/Cards/param_card.dat\n" # Change events in MG run card	
			# Generate MG events
			batchn += "./run.sh " + str(nevents) + " " + str(seeds[ij - ijobs[ibkg]]) +"\n" 
			
			# Munge LHE file
			lhe_file = "/tmp/sid_tmp/gridpack_" + str(mass) + "_" + str(ij) + "/events.lhe"
			batchn += "gunzip " + lhe_file + ".gz\n" #
			batchn += "cp " + lhe_file + " " + mg_dir + prefix + "/Events/events_" + str(mass) + "_" + str(ij) + ".lhe" "\n" #
			lhe_file = mg_dir + prefix + "/Events/events_" + str(mass) + "_" + str(ij) + ".lhe"

			batchn += "sed -i 's/49001010/4900101/g' "+lhe_file+"\n" # Change events in MG run card	
			batchn += "sed -i 's/49001011/4900101/g' "+lhe_file+"\n" # Change events in MG run card	
			batchn += "sed -i 's/49001012/4900101/g' "+lhe_file+"\n" # Change events in MG run card	
			batchn += "sed -i 's/49001013/4900101/g' "+lhe_file+"\n" # Change events in MG run card	
			batchn += "sed -i 's/49001014/4900101/g' "+lhe_file+"\n" # Change events in MG run card	

			# Clean up
			batchn += "rm -rf /tmp/sid_tmp/gridpack_" + str(mass) + "_" +  str(ij) + "\n" 

			# Run pythia for various rinv
			batchn += "cd " + base_dir + "/gen/\n" # Go to Pythia script folder
			for rinv in rinvs:
				batchn += "./monojet.exe -m lhe -w -i " + lhe_file + " -metmin 0 " + recluster_str + " -n " + str(nevents) + " -o " + data_dir + prefix + "/events_" + str(mass) + "_" + str(rinv) + "_" +  str(ij) +  " -nmatch " + str(nmatch[ibkg]) +" -inv " +str(rinv) + " -phimass 10 -lambda 5 -v \n"
			
			fname = "batch/batch_" + prefix + "_" + str(mass) + "_" + str(ij) + ".batch" # 
			f=open(fname, "w")
			f.write(batchn)
			f.close()
			os.system("sbatch " + fname);