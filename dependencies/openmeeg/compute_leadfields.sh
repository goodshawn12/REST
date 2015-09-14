#!/usr/bin/env bash

GEOMETRY=head_model.geom
CONDUCTIVITIES=head_model.cond
DIPOLES=cortex_dipoles.txt
EEG_ELECTRODES=eeg_channels_locations.txt
#SQUIDS=meg_channels_locations.squids
#INTERNAL_ELECTRODES=internal_electrodes_locations.txt

# Leadfields
EEG_LEADFIELD=eeg_leadfield.mat
#MEG_LEADFIELD=meg_leadfield.mat
#EIT_LEADFIELD=eit_leadfield.mat
#IP_LEADFIELD=ip_leadfield.mat

# Name temporary matrices
HM=tmp/hm.bin           # For EEG and MEG
HMINV=tmp/hm_inv.bin    # For EEG and MEG
DSM=tmp/dsm.bin         # For EEG and MEG
H2EM=tmp/h2em.bin       # For EEG
#H2MM=tmp/tmp.h2mm       # For MEG
#DS2MEG=tmp/tmp.ds2mm    # For MEG
#EITSM=tmp/tmp.eitsm     # For EIT
#IPHM=tmp/tmp.iphm       # For IP
#IPSM=tmp/tmp.ipsm       # For IP

mkdir -p tmp

# Compute EEG gain matrix
#for i in *.tri; do om_mesh_convert -i $i -o ${i/.tri/.tri} -invert ; done
#om_check_geom -g ${GEOMETRY}
om_assemble -HM ${GEOMETRY} ${CONDUCTIVITIES} ${HM}
om_minverser ${HM} ${HMINV}
om_assemble -DSM ${GEOMETRY} ${CONDUCTIVITIES} ${DIPOLES} ${DSM}
om_assemble -H2EM ${GEOMETRY} ${CONDUCTIVITIES} ${EEG_ELECTRODES} ${H2EM}
om_gain -EEG ${HMINV} ${DSM} ${H2EM} ${EEG_LEADFIELD}

# Compute MEG gain matrix
#om_assemble -H2MM ${GEOMETRY} ${CONDUCTIVITIES} ${SQUIDS} ${H2MM}
#om_assemble -DS2MM ${DIPOLES} ${SQUIDS} ${DS2MEG}
#om_gain -MEG ${HMINV} ${DSM} ${H2MM} ${DS2MEG} ${MEG_LEADFIELD}

# Compute EIT gain matrix
#om_assemble -EITSM ${GEOMETRY} ${CONDUCTIVITIES} ${EEG_ELECTRODES} ${EITSM}
#om_gain -EEG ${HMINV} ${EITSM} ${H2EM} ${EIT_LEADFIELD}

# Compute Internal Potential
#om_assemble -H2IPM ${GEOMETRY} ${CONDUCTIVITIES} ${INTERNAL_ELECTRODES} ${IPHM}
#om_assemble -DS2IPM ${GEOMETRY} ${CONDUCTIVITIES} ${DIPOLES} ${INTERNAL_ELECTRODES} ${IPSM}
#om_gain -IP ${HMINV} ${DSM} ${IPHM} ${IPSM} ${IP_LEADFIELD}

