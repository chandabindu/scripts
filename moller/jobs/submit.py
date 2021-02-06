import sys
import os
import subprocess
import time

if not sys.argv[1]:
        sourceDir = "/lustre19/expphy/volatile/halla/moller12gev/chandan/sim_out"
        #sourceDir = "/lustre19/expphy/volatile/halla/parity/chandan/sim_out"
else:
        sourceDir = sys.argv[1]
if not sys.argv[2]:
        config = "AsymmetricCoil4_ep"
else:
        config = sys.argv[2]

if not sys.argv[3]:
        generator = "ep"
else:
        generator = sys.argv[3]

runrange= range(int(sys.argv[4]),int(sys.argv[5]))

jobs=sourceDir+"/"+config+"/jobs"
outDir=sourceDir+"/"+config+"/output"
if not os.path.exists(jobs):
  os.system("mkdir "+jobs)
if not os.path.exists(outDir):
  os.system("mkdir "+outDir)


home = sourceDir+"/"+config

for i in runrange:
  filename=sourceDir+"/"+config+"/remollout_"+generator+str(i)+".root"
  outfile=outDir+"/analysis_"+str(i)+".root"
  if os.path.exists(filename):
    jsubf=open(jobs+"/ep_"+str(i)+".sh", "w")
    jsubf.write("#!/bin/bash\n")
    jsubf.write("#SBATCH --account=halla\n")
    jsubf.write("#SBATCH --partition=production\n")
    jsubf.write("#SBATCH --job-name=remollAna\n")
    jsubf.write("#SBATCH --time=00:05:00\n")
    jsubf.write("#SBATCH --nodes=1\n")
    jsubf.write("#SBATCH --ntasks=1\n")
    jsubf.write("#SBATCH --cpus-per-task=1\n")
    jsubf.write("#SBATCH --mem=3G\n")
    jsubf.write("tcsh -c \"source /apps/root/6.18.00/setroot_CUE\"\n")
    jsubf.write("#SBATCH --chdir="+outDir+"\n")
    #jsubf.write("#SBATCH --output="+outDir+"/tmp"+"_"+str(i)+".out\n")
    #jsubf.write("#SBATCH --output=/lustre19/expphy/volatile/halla/moller12gev/chandan/sim_out/OffsetCoil3_ep/jobs/out.out\n")
    jsubf.write("#SBATCH --output=/lustre19/expphy/volatile/halla/parity/chandan/sim_out/OffsetCoil3_ep/jobs/out.out\n")
    #jsubf.write("#SBATCH --error="+outDir+"/tmp"+"_"+str(i)+".err\n")
    jsubf.write("cd "+home+"\n")
    jsubf.write("echo \"Current working directory is `pwd`\"\n")	
    jsubf.write("./reroot -q -b ShiftedCoil_BoreVsAcceptance.C'(\""+filename+"\", \""+outfile+"\")'\n")
    print("sbatch "+jobs+"/ep_"+str(i)+".sh") 
      		
    
		
	
	
