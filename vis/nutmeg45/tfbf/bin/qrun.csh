#!/bin/csh
#
#   1.  Define resource requirements:
#       -q   =>queue_name
#       -t   => number of tasks
#       -cwd => run from cwd
#$ -q bil.q
#$ -l arch=lx24-x86
#$ -t 1
#$ -cwd
#$ -j y
#$ -o qlogs
#$ -e qlogs
#$ -r y

hostname
source /data/research_meg/tfbf/bin/MatlabSetup.csh
time $*
