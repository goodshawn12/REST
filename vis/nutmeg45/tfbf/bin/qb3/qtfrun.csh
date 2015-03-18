#!/bin/csh
#   1.  Define resource requirements:
#       -q   =>queue_name
#       -t   => number of tasks
#       -cwd => run from cwd
#$ -cwd
#$ -l intel64=false
##$ -l opt64=true
#$ -l h_rt=00:40:00
#$ -l mem_free=600M
#$ -j y
#$ -o qlogs
#$ -e qlogs
#$ -r y
#$ -q medium.q

hostname
source ~/bin/MatlabSetup.csh

time ~/bin/`arch`/tfrun $1 $2 $SGE_TASK_ID $3

