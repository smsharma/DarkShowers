import sys, os
import numpy as np
from math import log10, floor

batch='''#!/bin/bash
#SBATCH -N 1   # node count
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4gb
#SBATCH -t 20:00:00
##SBATCH --mail-type=begin
##SBATCH --mail-type=end
##SBATCH --mail-user=smsharma@princeton.edu
#SBATCH -p hepheno

cd /group/hepheno/smsharma/Dark-Showers
source env.sh
cd /group/hepheno/smsharma/Dark-Showers/gen
source activate venv_py27

./monojet.exe -m lhe -w -i '''

# rinv_ary = np.linspace(0.1,1.,5)
rinv_ary = [0.01,0.1,.2,.3,.4,.5,.6,.7,.8,.9, .98, 0.99, 1.0]

# for MZp in [500, 1000, 3000]:
for MZp in range(400,3000,200):
    for rinv in rinv_ary:    
        fname_tag = "/group/hepheno/hlou/resonance_" + str(MZp) + ".lhe -metmin 0 -n 100000 -phimass 10 -lambda 5 -v "
        rinv_tag = "-inv " + str(rinv)
        out_tag = "-o " + "Resonance/Resonance_M"+str(MZp)+"_rinv" + str(rinv)
        tag = fname_tag + " " + rinv_tag + " " + out_tag
        print tag
        batchn = batch + tag
        fname = "batch/run_monojet_"+str(MZp)+"_"+str(rinv)+".batch"
        f=open(fname, "w")
        f.write(batchn)
        f.close()
        os.system("sbatch "+fname);
