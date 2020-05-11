"""
Code to generate large amount of parton-level
information in MadGraph and pipe it through a
Pythia script.
"""

import sys, os
import random

import numpy as np

# extra_tag = 'test1new4flxqptjptj1_'
extra_tag = ''
base_dir = "/group/hepheno/smsharma/Dark-Showers/"
# mg_base = "MG5_aMC_v2_5_6_patch_signal"
mg_base = "MG5_aMC_v2_5_6_patch_int_signal"
# mg_base = "MG5_aMC_v2_5_5"
mg_dir = base_dir + mg_base # MadGraph directory
data_dir = base_dir + "/data/"

# Which signals should generate (gridpacks should exist in these dirs already)
bkg_prefix = ["sig_tchannel", "sig_zprime"]

mdark = 10 # Dark quark mass in GeV

nevents = 100000 # Events per job
ijobs = [0, 0] # Start index
njobs = [1, 0] # Number of jobs to submit
 # Prefix that the MG job was generated with
nmatch = [2, 2] # Pythia nmatch parameter
recluster = False # Whether to recluster jets in pythia
hepmc = False
# higher_masses = [3000,5000,10000,20000,50000,100000]
# higher_masses = [100000]
# masses = [list(np.arange(500,3000,100)) + higher_masses,list(np.arange(500,5000,100)) + higher_masses[1:]]
# masses = [list(np.arange(500,2500,100)),list(np.arange(500,5000,100)) + higher_masses[1:]]
# masses = [list(np.arange(500,3000,100)),list(np.arange(500,5000,100))]
# masses = [[],list(np.arange(500,5000,100))]
# masses = [list(np.arange(1000,3000,100)) + higher_masses,list(np.arange(500,3000,100))]
# masses = [list(np.arange(500,3000,100)) + higher_masses,[]]
# masses = [[],list(np.arange(500,5000,100)) + higher_masses[1:]]
# masses = [[],list(np.arange(500,5000,100)) + higher_masses[1:]]
# masses = [list(np.arange(500,3500,1000)),list(np.arange(500,3500,1000))]
# masses = [[],list(np.arange(3000,5100,100))]
# masses = [higher_masses,higher_masses]
# masses = [[20000,50000,100000],[]]
masses = [list(np.arange(100,2000, 100)) + [100000],[]]
# masses = [higher_masses,[]]
# masses = [[],[100000]]

rinvs = [0.0, 0.01,0.1,.2,.3,.4,.5,.6,.7,.8,.9, 0.99, 1.0]

if recluster:
	recluster_str = " -Zprime "
else:
	recluster_str = " "

batch='''#!/bin/bash
#SBATCH -N 1   # node count
#SBATCH --ntasks-per-node=6
#SBATCH --mem=24gb
#SBATCH -t 48:00:00
##SBATCH --mail-type=begin
##SBATCH --mail-type=end
##SBATCH --mail-user=smsharma@princeton.edu
#SBATCH -p dept
source /opt/rh/devtoolset-3/enable
source activate venv_py27
source /group/hepheno/smsharma/Dark-Showers/env.sh
source /group/hepheno/heptools/root/bin/thisroot.sh
'''

for ibkg, prefix in enumerate(bkg_prefix): 

	print "Submitting to generate", nevents*njobs[ibkg], "events for", prefix, "! Nice"

	# Make folder for final ntuples
	if not os.path.exists(data_dir + prefix):
		os.mkdir(data_dir + prefix)

	for imass, mass in enumerate(masses[ibkg]):

		seeds = random.sample(range(1, 30081), njobs[ibkg])

		for ij in range(ijobs[ibkg], ijobs[ibkg] + njobs[ibkg]):
			
			batchn = batch + mg_dir + prefix + "\n" # Go to appropriate MG folder
			
			# Copy over gridpack to temp
			batchn += "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/group/hepheno/heptools/HepMC-2.06.09/lib/\n" # Clean up previous gridpack
			# batchn += "rm -r /tmp/sid_tmp/MG5_aM*\n" # Clean up previous gridpack
			# batchn += "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/group/hepheno/heptools/HepMC-2.06.09/lib/\n" # Clean up previous gridpack
			# batchn += "rm -rf /tmp/sid_tmp/" + mg_base + "_"+ prefix + "_" + str(mass) + "_" +  str(ij) +"\n" # Clean up previous gridpack
			# batchn += "mkdir /tmp/sid_tmp/" +"\n" # Make folder for this run
			batchn += "cp -r "+ mg_dir + " " + base_dir + mg_base + "_"+ prefix + "_" + str(mass) + "_"  + str(ij) +"\n" # Copy over gridpack files
			batchn += "cd " + base_dir + mg_base + "_"+ prefix + "_" + str(mass) + "_"  + str(ij) + "/bin/" + prefix + "\n" # Go to folder
			batchn += "rm RunWeb" +"\n" # Delete random shit file if needed			
			batchn += "sed -i 's/<PhiMass>/"+str(mass)+"/g' Cards/param_card.dat\n" # Change events in MG run card	
			batchn += "sed -i 's/<DarkQuarkMass>/"+str(mdark)+"/g' Cards/param_card.dat\n" # Change events in MG run card	
			batchn += "sed -i 's/<Seed>/"+str(seeds[ij - ijobs[ibkg]])+"/g' Cards/run_card.dat\n" # Change events in MG run card	 
			# Generate MG events
			batchn += "./bin/generate_events -f " +"\n" 

			# Munge LHE file
			lhe_file = base_dir + mg_base + "_"+ prefix + "_" + str(mass) + "_"  + str(ij) + "/bin/"+prefix+"/Events/run_01/unweighted_events.lhe"
			batchn += "gunzip " + lhe_file + ".gz\n" #
			# batchn += "cp " + lhe_file + " " + mg_dir + "/bin/" + prefix + "/Events/events_"+extra_tag+  str(mass) + "_" + str(ij) + ".lhe" "\n" #
			
			# lhe_file = mg_dir + "/bin/" + prefix + "/Events/events_"+extra_tag+  str(mass) + "_" + str(ij) + ".lhe"

			# if prefix == "sig_tchannel":
			# 	batchn += "sed -i 's/49001010/4900101/g' "+lhe_file+"\n" # Change events in MG run card	
			# 	batchn += "sed -i 's/49001011/4900101/g' "+lhe_file+"\n" # Change events in MG run card	
			# 	batchn += "sed -i 's/49001012/4900101/g' "+lhe_file+"\n" # Change events in MG run card	
			# 	batchn += "sed -i 's/49001013/4900101/g' "+lhe_file+"\n" # Change events in MG run card	
			# 	batchn += "sed -i 's/49001014/4900101/g' "+lhe_file+"\n" # Change events in MG run card	
			# elif prefix == "sig_zprime":
			# 	batchn += "sed -i 's/5000521/4900101/g' "+lhe_file+"\n" # Change events in MG run card	

			# # Clean up
			# batchn += "rm -rf /tmp/sid_tmp/" + mg_base + "_"+ prefix + "_" + str(mass) + "_" +  str(ij) +"\n" # Clean up previous gridpack

			# # Run pythia for various rinv
			# batchn += "cd " + base_dir + "/gen/\n" # Go to Pythia script folder
			
			# if hepmc:
			# 	batchn += 'mkdir ' + data_dir + prefix + '/hepmc \n'
			# 	hepmc_str = ' -hepmc ' + data_dir + prefix + '/hepmc/events_hepmc_M' + str(mass) + '_rinv\${rinv}_' +  str(ij)
			# else:
			# 	hepmc_str = ' '

			# batchn += 'echo "#!/bin/bash \nrinvs=(0.0 0.01 0.1 .2 .3 .4 .5 .6 .7 .8 .9 0.99 1.0)\necho i = \$1 \nexport rinv=\${rinvs[\$1]} \n./monojet.exe -m lhe -w -i ' + lhe_file + ' -metmin 0 ' + recluster_str + ' -n ' + str(nevents) + hepmc_str +  ' -o ' + data_dir + prefix + '/events_' +extra_tag+ str(mass) + '_\${rinv}_' +  str(ij) +  ' -nmatch ' + str(nmatch[ibkg]) + ' -inv \$rinv -phimass ' + str(mdark) + ' -lambda ' + str(mdark/2.) + ' -v " > run-'+prefix+'_' + extra_tag+str(mass)+'_'+str(ij)+'.sh\n'
			# batchn += 'chmod u+x run-'+prefix+'_' + extra_tag+str(mass)+'_'+str(ij)+'.sh\n'

			# batchn += 'echo  0-12 ./run-'+prefix+'_' + extra_tag + str(mass)+'_'+str(ij)+'.sh %t  > run-'+prefix+'_' + extra_tag+str(mass)+'_'+str(ij)+'.conf\n'+'\n'+'\n'+'srun --multi-prog --no-kill --wait=0 run-'+prefix+'_' + extra_tag+str(mass)+'_'+str(ij)+'.conf'+'\n'+'\n'
			
			fname = "batch/batch_" + prefix + "_" + extra_tag + str(mass) + "_" + str(ij) + ".batch" # 
			f=open(fname, "w")
			f.write(batchn)
			f.close()
			os.system("sbatch " + fname);