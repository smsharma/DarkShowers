"""
Code to generate large amount of parton-level
information in MadGraph.
"""

import sys, os
import random

import numpy as np

extra_tag = "md10_"
base_dir = "/group/hepheno/smsharma/Dark-Showers/"
mg_base = "MG5_aMC_v2_5_6_patch_int_signal"
mg_dir = base_dir + mg_base  # MadGraph directory

# Which signals should generate (gridpacks should exist in these dirs already)
sig_prefix = ["sig_tchannel_lambda1", "sig_tchannel_lambda0p3"]

mdark = 10  # Dark quark mass in GeV

nevents = 10000  # Events per job
ijobs = [0, 0]  # Start indices of job
njobs = [16, 16]  # Number of jobs to submit for each signal
masses = [list(np.arange(100, 2000, 100)) + [100000] + list(np.arange(2000, 10000, 500))]  # Mediator masses to loops over

batch = """#!/bin/bash
#SBATCH -N 1   # node count
#SBATCH --ntasks-per-node=13
#SBATCH --mem=47gb
#SBATCH -t 3:00:00
##SBATCH --mail-type=begin
##SBATCH --mail-type=end
##SBATCH --mail-user=smsharma@princeton.edu
#SBATCH -p hepheno
source /group/hepheno/smsharma/Dark-Showers/env.sh
source /group/hepheno/heptools/root/bin/thisroot.sh
"""

for isig, prefix in enumerate(sig_prefix):

    print "Submitting to generate", nevents * njobs[isig], "events for", prefix, "! Nice"

    for imass, mass in enumerate(masses[isig]):

        print "Submitting", mass

        seeds = random.sample(range(1, 30081), njobs[isig])

        for ij in range(ijobs[isig], ijobs[isig] + njobs[isig]):

            batchn = batch + mg_dir + prefix + "\n"  # Go to appropriate MG folder

            ## Copy over gridpack to temp
            batchn += "rm -rf /tmp/sid_tmp/" + mg_base + "_" + prefix + "_" + str(mass) + "_" + str(ij) + "\n"  # Clean up previous gridpack
            batchn += "mkdir /tmp/sid_tmp/" + "\n"  # Make folder for this run
            batchn += "cp -r " + mg_dir + " /tmp/sid_tmp/" + mg_base + "_" + prefix + "_" + str(mass) + "_" + str(ij) + "\n"  # Copy over gridpack files
            batchn += "cd /tmp/sid_tmp/" + mg_base + "_" + prefix + "_" + str(mass) + "_" + str(ij) + "/bin/" + prefix + "\n"  # Go to folder
            batchn += "rm RunWeb" + "\n"

            batchn += "sed -i 's/<PhiMass>/" + str(mass) + "/g' Cards/param_card.dat\n"  # Change params in MG card
            batchn += "sed -i 's/<DarkQuarkMass>/" + str(mdark) + "/g' Cards/param_card.dat\n"  # Change params in MG card
            batchn += "sed -i 's/<Seed>/" + str(seeds[ij - ijobs[isig]]) + "/g' Cards/run_card.dat\n"  # Change params in MG card

            # Generate MG events
            batchn += "./bin/generate_events -f " + "\n"

            # Munge LHE file
            lhe_file = "/tmp/sid_tmp/" + mg_base + "_" + prefix + "_" + str(mass) + "_" + str(ij) + "/bin/" + prefix + "/Events/run_01/unweighted_events.lhe"
            html_file = "/tmp/sid_tmp/" + mg_base + "_" + prefix + "_" + str(mass) + "_" + str(ij) + "/bin/" + prefix + "/HTML/run_01/results.html"
            batchn += "gunzip " + lhe_file + ".gz\n"  #
            batchn += "cp " + lhe_file + " " + mg_dir + "/bin/" + prefix + "/Events/events_" + extra_tag + str(mass) + "_" + str(ij) + ".lhe" "\n"  #
            batchn += "cp " + html_file + " " + mg_dir + "/bin/" + prefix + "/Events/html_" + extra_tag + str(mass) + "_" + str(ij) + ".html" "\n"  #

            lhe_file = mg_dir + "/bin/" + prefix + "/Events/events_" + extra_tag + str(mass) + "_" + str(ij) + ".lhe"

            if prefix in ["sig_tchannel_lambda1", "sig_tchannel_lambda0p3"]:
                batchn += "sed -i 's/49001010/4900101/g' " + lhe_file + "\n"  # Change events in MG run card
                batchn += "sed -i 's/49001011/4900101/g' " + lhe_file + "\n"  # Change events in MG run card
                batchn += "sed -i 's/49001012/4900101/g' " + lhe_file + "\n"  # Change events in MG run card
                batchn += "sed -i 's/49001013/4900101/g' " + lhe_file + "\n"  # Change events in MG run card
                batchn += "sed -i 's/49001014/4900101/g' " + lhe_file + "\n"  # Change events in MG run card
            elif prefix == "sig_zprime":
                batchn += "sed -i 's/5000521/4900101/g' " + lhe_file + "\n"  # Change events in MG run card

            # Clean up
            batchn += "rm -rf /tmp/sid_tmp/" + mg_base + "_" + prefix + "_" + str(mass) + "_" + str(ij) + "\n"  # Clean up previous gridpack

            fname = "batch/batch_" + prefix + "_" + extra_tag + str(mass) + "_" + str(ij) + ".batch"  #
            f = open(fname, "w")
            f.write(batchn)
            f.close()
            os.system("sbatch " + fname)

