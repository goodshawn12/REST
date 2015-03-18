#!/bin/csh
#   1.  Define resource requirements:
#       -q   =>queue_name
#       -t   => number of tasks
#       -cwd => run from cwd
#$ -cwd
#$ -l intel64=false
##$ -l h_rt=00:4:00
#$ -l mem_free=1200M
#$ -j y
#$ -o qlogs
#$ -e qlogs
#$ -r y
#$ -q short.q

hostname
source ~/bin/MatlabSetup.csh

~/bin/`arch`/nut_tfbf2timef $*
