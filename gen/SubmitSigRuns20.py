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
#SBATCH -p hepheno

cd /group/hepheno/smsharma/Dark-Showers
source env.sh
cd /group/hepheno/smsharma/Dark-Showers/gen
source activate venv_py27

./monojet.exe -m lhe -w -i /group/hepheno/heptools/MG5_aMC_v2_5_2/bin/DMsimp_s1_20/Events/run_01/unweighted_events.lhe -metmin 0 -n 500000 -phimass 20 -lambda 10 '''

# rinv_ary = np.linspace(0.1,1.,5)
rinv_ary = [0.01,0.1,.2,.3,.4,.5,.6,.7,.8,.9, .98, 0.99, 1.0]

for rinv in rinv_ary:    
    rinv_tag = "-inv " + str(rinv)
    out_tag = "-o " + "monojet_20_rinv" + str(rinv)
    tag = rinv_tag + " " + out_tag
    print tag
    batchn = batch + tag
    fname = "batch/run_monojet_20_"+str(rinv)+".batch"
    f=open(fname, "w")
    f.write(batchn)
    f.close()
    os.system("sbatch "+fname);