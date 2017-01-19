import sys, os
import numpy as np
from math import log10, floor

rtn = lambda x, n: round(x, -int(floor(log10(x))) + (n - 1))

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

./monojet.exe -n 10000'''

rinv_ary = np.linspace(0.1,1.,5)
mbf_ary = np.logspace(np.log10(100),np.log10(5000),5)

for rinv in rinv_ary:    
    for mbf in mbf_ary:
        mbf_tag = "-mphi " + str(mbf)
        rinv_tag = "-inv " + str(rinv)
        out_tag = "-o " + "tch_nomet_mphi" + str(mbf) + "_rinv" + str(rinv)
        tag = mbf_tag + " " + rinv_tag + " " + out_tag
        print tag
        batchn = batch + tag
        fname = "batch/run_tch_nomet_"+str(mbf)+"_"+str(rinv)+".batch"
        f=open(fname, "w")
        f.write(batchn)
        f.close()
        os.system("sbatch "+fname);