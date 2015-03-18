#!/bin/csh
#$ -q bil.q
#$ -l arch=lx24-x86
#$ -l mem_free=1200M
#$ -cwd
#$ -j y
#$ -o qlogs
#$ -e qlogs
#$ -r y
#$ -V

# 3/13/07 AMF modified to remove paths from the script.
# Thus, you should have copies of MatlabSetup.csh (from
# /data/research_meg/tfbf/bin/MatlabSetup.csh) in your
# data directory, along with your parameter files (300ms.mat,
# etc.

# '$argv[2-]' takes the last two inputs from qTFBF.

# Also, if you have very high memory intensive data,
# please increase mem_free to something like 1800M.

# tfcov is in /data/research_meg/tfbf/bin/i686, 
# which should be in your path,
# or /data/research-meg/tfbf/bin/`arch`/tfcov

hostname
source MatlabSetup.csh
# from /data/research_meg/tfbf/bin/MatlabSetup.csh


#4-12 Hz
if ($SGE_TASK_ID == 1) then
    time tfcov $1 300ms.mat 4 12 $argv[2-]
endif

#12-30 Hz
if ($SGE_TASK_ID == 2) then
    time tfcov $1 200ms.mat 12 30 $argv[2-]
endif

#30-55 Hz
if ($SGE_TASK_ID == 3) then
    time tfcov $1 150ms.mat 30 55 $argv[2-]
endif

#65-90 Hz
if ($SGE_TASK_ID == 4) then
    time tfcov $1 100ms.mat 65 90 $argv[2-]
endif

#90-115 Hz
if ($SGE_TASK_ID == 5) then
    time tfcov $1 100ms.mat 90 115 $argv[2-]
endif

#125-150 Hz
if ($SGE_TASK_ID == 6) then
    time tfcov $1 100ms.mat 125 150 $argv[2-]
endif

#150-175 Hz
if ($SGE_TASK_ID == 7) then
    time tfcov $1 100ms.mat 150 175 $argv[2-]
endif

#185-300 Hz
if ($SGE_TASK_ID == 8) then
    time tfcov $1 100ms.mat 185 300 $argv[2-]
endif
