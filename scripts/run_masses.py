import os
# import paramiko
import numpy as np

batch='''#!/bin/sh
#PBS -j oe
#PBS -o /home/smsharma/scratch/job.out
#PBS -l nodes=1:ppn=1 -l cput=100:59:00
#PBS -q dept
#PBS -m abe
#PBS -k oe
#PBS -M smsharma@princeton.edu

cd /group/hepheno/venv
source bin/activate
source /group/hepheno/smsharma/Dark_Showers/env.sh
cd /group/hepheno/smsharma/Dark_Showers/gen
source /group/hepheno/heptools/root/bin/thisroot.sh
./jets_simple.exe -m tchannel -n 5000 -mphi '''

mDMList = np.array([100,400,600,800,1000,2000,3000])
for mDM in mDMList:
    batchn = batch + str(mDM) + ' -o tchannel_mchitilde_' + str(mDM)
    fname = "run_mphi_"+str(mDM)+".batch"
    f=open(fname, "w")
    f.write(batchn)
    f.close()
    os.system("qsub -l mem=16gb "+fname);
