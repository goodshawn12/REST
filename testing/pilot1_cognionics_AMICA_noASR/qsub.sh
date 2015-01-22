#!/bin/bash
#$ -cwd
#$ -q q5
#$ -j y
#$ -S /bin/bash
#$ -pe mpich 24
/home/jason/mpich2-1.5-install/bin/mpirun -np $NSLOTS -machinefile $TMP/machines /data/common/matlab/eeglab/plugins/amica1.0/amica15_c /home/lpiontonachini/Dropbox/School/Research/VisEEG_local/AMICA_noASR/input.param