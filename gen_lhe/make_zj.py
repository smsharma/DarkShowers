import sys, os
import random

n_jobs = 4
prefix = "bkg_Zj"
nevt_each = 500000

batch='''#!/bin/bash
#SBATCH -N 1   # node count
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16gb
#SBATCH -t 16:00:00
##SBATCH --mail-type=begin
##SBATCH --mail-type=end
##SBATCH --mail-user=smsharma@princeton.edu
#SBATCH -p hepheno
source activate venv_py27
cd /group/hepheno/heptools/MG5_aMC_v2_5_2/bin/'''

for i in range(n_jobs):
	os.system("cp -r /group/hepheno/heptools/MG5_aMC_v2_5_2/bin/" 
		+ str(prefix) + 
		" /group/hepheno/heptools/MG5_aMC_v2_5_2/bin/"
		+ prefix+"_" + str(i))
	batchn = batch + prefix + "_" + str(i) + "/bin/\n"
	seed = random.randrange(0,123120) # Seed
	batchn += "sed -i 's/0   = iseed/"+str(seed)+"   = iseed/g' ../Cards/run_card.dat\n"
	batchn += "sed -i 's/50000 = nevents/"+str(nevt_each)+" = nevents/g' ../Cards/run_card.dat\n"
	batchn += "rm -r /group/hepheno/heptools/MG5_aMC_v2_5_2/bin/"+prefix+"_" + str(i)+"/RunWeb\n"
	batchn += "./generate_events -f\n"
	batchn += "gunzip /group/hepheno/heptools/MG5_aMC_v2_5_2/bin/" + prefix+"_" + str(i) +"/Events/run_01/unweighted_events.lhe.gz\n"
	batchn += "cd /group/hepheno/smsharma/Dark-Showers/gen/\n"
	batchn += "./monojet.exe -m lhe -w -i /group/hepheno/heptools/MG5_aMC_v2_5_2/bin/" + prefix + "_" + str(i) + "/Events/run_01/unweighted_events.lhe -metmin 0 -n 2000000 -o bkg_Zj_" + str(i) + " -nmatch 2\n"
	fname = "batch/batch_" + prefix + "_" + str(i) + ".batch"
	f=open(fname, "w")
	f.write(batchn)
	f.close()
	os.system("sbatch " + fname);
