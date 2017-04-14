import sys, os
import numpy as np
from math import log10, floor

batch='''#!/bin/bash
#SBATCH -N 1   # node count
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16gb
#SBATCH -t 20:00:00
##SBATCH --mail-type=begin
##SBATCH --mail-type=end
##SBATCH --mail-user=smsharma@princeton.edu
##SBATCH -p hepheno

cd /group/hepheno/smsharma/Dark-Showers
source env.sh
cd /group/hepheno/smsharma/Dark-Showers/gen
source activate venv_py27

./monojet.exe -m lhe -w -i '''

# rinv_ary = np.linspace(0.1,1.,5)
rinv_ary = [0.01,0.1,.2,.3,.4,.5,.6,.7,.8,.9, .98, 0.99, 1.0]

for Mq in [100, 200, 500, 1000]:
    for rinv in rinv_ary:    
        fname_tag = "/group/hepheno/heptools/MG5_aMC_v2_5_2/bin/DMsimp_s1_" + str(Mq) + "/Events/run_01/unweighted_events.lhe -metmin 0 -n 200000 -phimass " +str(2*Mq) + " -lambda " + str(Mq) + " -v "
        rinv_tag = "-inv " + str(rinv)
        out_tag = "-o " + "monojet_"+str(Mq)+"_rinv" + str(rinv)
        tag = fname_tag + " " + rinv_tag + " " + out_tag
        print tag
        batchn = batch + tag
        fname = "batch/run_monojet_"+str(Mq)+"_"+str(rinv)+".batch"
        f=open(fname, "w")
        f.write(batchn)
        f.close()
        os.system("sbatch "+fname);
