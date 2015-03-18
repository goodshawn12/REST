#!/bin/csh
# -q bil.q

#$ -q bil.q@bodha.radiology.ucsf.edu
#$ -q bil.q@advaita.radiology.ucsf.edu
#$ -q bil.q@carvaka.radiology.ucsf.edu
#$ -q bil.q@kaumudi.radiology.ucsf.edu
#$ -q bil.q@bil2.radiology.ucsf.edu
#$ -q bil.q@vedanta.radiology.ucsf.edu
#$ -q bil.q@nirukta.radiology.ucsf.edu
# -q bil.q@sidhanta.radiology.ucsf.edu
#$ -q bil.q@tarka.radiology.ucsf.edu
#$ -q bil.q@ahimsa.radiology.ucsf.edu
# -q bil.q@mimamsa.radiology.ucsf.edu




#$ -cwd
#$ -l mem_free=1200M
#$ -l arch=lx24-x86
#$ -j y
#$ -o qlogs
#$ -e qlogs
#$ -r y
#$ -V

# setenv SGE_TASK_ID covwindownum  % set SGE_TASK_ID if not running on SGE cluster
# ./qtfrun.csh sessionfile outputfile algo


hostname
source /data/research_meg/tfbf/bin/MatlabSetup.csh

time /data/research_meg/tfbf/bin/`arch`/tfrun $1 $2 $SGE_TASK_ID $3

