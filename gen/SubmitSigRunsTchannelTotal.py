import sys, os
import numpy as np
from math import log10, floor

batch='''#!/bin/bash
#SBATCH -N 1   # node count
#SBATCH --ntasks-per-node=1
#SBATCH --mem=3gb
#SBATCH -t 2:00:00
##SBATCH --mail-type=begin
##SBATCH --mail-type=end
##SBATCH --mail-user=smsharma@princeton.edu
#SBATCH -p hepheno

cd /group/hepheno/smsharma/Dark-Showers
source env.sh
cd /group/hepheno/smsharma/Dark-Showers/gen
source activate venv_py27

./monojet.exe -m lhe -w -Zprime -i '''

# rinv_ary = np.linspace(0.1,1.,5)
rinv_ary = [0.01,0.1,.2,.3,.4,.5,.6,.7,.8,.9, .98, 0.99, 1.0]

# for MZp in [500, 1000, 3000]:
# for MZp in [500,1000,5000,8000]:
# masses = np.arange(2500,4500,500)
# for i in :
# masses = [500,1000,2000,3000,4000,5000,7000,10000,20000,50000,100000]
masses = [100000]
masses = [500,1000,2000,3000,4000]
for MZp in masses:
    for rinv in rinv_ary:    
        fname_tag = "/group/hepheno/smsharma/Dark-Showers/MG5_aMC_v2_5_2/bin/tchannel_total_test_" + str(MZp) + "/Events/run_01/unweighted_events.lhe -metmin 0 -n 200000 -phimass 10 -lambda 5 -v "
        rinv_tag = "-inv " + str(rinv)
        out_tag = "-o " + "tChannelDirect/tchannel_total_M"+str(MZp)+"_rinv" + str(rinv)
        tag = fname_tag + " " + rinv_tag + " " + out_tag
        print tag
        batchn = batch + tag
        fname = "batch/run_monojet_"+str(MZp)+"_"+str(rinv)+".batch"
        f=open(fname, "w")
        f.write(batchn)
        f.close()
        os.system("sbatch "+fname);
