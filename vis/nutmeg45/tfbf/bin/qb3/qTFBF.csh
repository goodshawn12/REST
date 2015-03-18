#!/bin/csh
#   1.  Define resource requirements:
#       -q   =>queue_name
#       -t   => number of tasks
#       -cwd => run from cwd
##$ -t 1-8
#$ -cwd
##$ -l opt64=true
#$ -l intel64=false
##$ -l h_rt=00:25:00
#$ -l mem_free=1200M
#$ -j y
#$ -o qlogs
#$ -e qlogs
#$ -r y
#$ -q short.q

hostname
source ~/bin/MatlabSetup.csh

#4-12 Hz
if ($SGE_TASK_ID == 1) then
    time ~/bin/`arch`/tfcov $1 300ms.mat 4 12 $argv[2-]
endif

#12-30 Hz
if ($SGE_TASK_ID == 2) then
    time ~/bin/`arch`/tfcov $1 200ms.mat 12 30 $argv[2-]
endif

#30-55 Hz
if ($SGE_TASK_ID == 3) then
    time ~/bin/`arch`/tfcov $1 150ms.mat 30 55 $argv[2-]
endif

#65-90 Hz
if ($SGE_TASK_ID == 4) then
    time ~/bin/`arch`/tfcov $1 100ms.mat 65 90 $argv[2-]
endif

#90-115 Hz
if ($SGE_TASK_ID == 5) then
    time ~/bin/`arch`/tfcov $1 100ms.mat 90 115 $argv[2-]
endif

#125-150 Hz
if ($SGE_TASK_ID == 6) then
    time ~/bin/`arch`/tfcov $1 100ms.mat 125 150 $argv[2-]
endif

#150-175 Hz
if ($SGE_TASK_ID == 7) then
    time ~/bin/`arch`/tfcov $1 100ms.mat 150 175 $argv[2-]
endif

#185-300 Hz
if ($SGE_TASK_ID == 8) then
    time ~/bin/`arch`/tfcov $1 100ms.mat 185 300 $argv[2-]
endif
