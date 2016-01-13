import os
# import paramiko
import numpy as np

batch='''#!/bin/sh
#PBS -j oe
#PBS -o /home/smsharma/scratch/job.out
#PBS -l nodes=1:ppn=1 -l cput=20:59:00
#PBS -q hepheno
#PBS -m abe
#PBS -k oe
#PBS -M smsharma@princeton.edu

cd /group/hepheno/venv
source bin/activate
source /group/hepheno/smsharma/Dark_Showers/env.sh
cd /group/hepheno/smsharma/Dark_Showers/gen
source /group/hepheno/heptools/root/bin/thisroot.sh
./jets_simple.exe -m tchannel -n 1000 -alpha 0.1 -ptcut 600 -run 0 -mphi '''

mDMList = np.array([100,400,600,800,1000,2000,3000])
for mDM in mDMList:
    batchn = batch + str(mDM) + ' -o tchannel_mchitilde_alpha_' + str(mDM)
    fname = "run_mphi_"+str(mDM)+".batch"
    f=open(fname, "w")
    f.write(batchn)
    f.close()
    os.system("qsub -l mem=4gb "+fname);
