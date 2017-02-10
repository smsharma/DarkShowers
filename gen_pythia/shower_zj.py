import sys, os
import random

n_jobs = 10
prefix = "bkg_Zj"

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
cd /group/hepheno/smsharma/Dark-Showers/gen/
./monojet.exe -m "lhe" -w -i /group/hepheno/heptools/MG5_aMC_v2_5_2/bin/'''

for i in range(n_jobs):
	batchn = batch + prefix + "/Events/run_01/unweighted_events.lhe -metmin 0 -n 2000000 -o bkg_Zj -nmatch 2\n"
	fname = "batch/batch_" + prefix + "_" + str(i) + ".batch"
	f=open(fname, "w")
	f.write(batchn)
	f.close()
	os.system("sbatch " + fname);