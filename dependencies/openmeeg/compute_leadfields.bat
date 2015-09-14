set GEOMETRY=head_model.geom
set CONDUCTIVITIES=head_model.cond
set DIPOLES=cortex_dipoles_small.txt
set EEG_ELECTRODES=eeg_channels_locations.txt
set SQUIDS=meg_channels_locations.squids
set INTERNAL_ELECTRODES=internal_electrodes_locations.txt

rem Leadfields
set EEG_LEADFIELD=eeg_leadfield.mat
set MEG_LEADFIELD=meg_leadfield.mat
set EIT_LEADFIELD=eit_leadfield.mat
set IP_LEADFIELD=ip_leadfield.mat

rem Name temporary matrices
rem For EEG and MEG

set HM=tmp\tmp.hm
set HMINV=tmp\tmp.hm_inv
set DSM=tmp\tmp.dsm

rem For EEG
set H2EM=tmp\tmp.h2em

rem for MEG
set H2MM=tmp\tmp.h2mm
set DS2MEG=tmp\tmp.ds2mm

rem for EIT
set EITSM=tmp\tmp.eitsm

rem for IP (Internal Potential)
set IPHM=tmp\tmp.iphm
set IPSM=tmp\tmp.ipsm

mkdir tmp

rem Compute EEG gain matrix
om_assemble -HM %GEOMETRY% %CONDUCTIVITIES% %HM%
om_minverser %HM% %HMINV%
om_assemble -DSM %GEOMETRY% %CONDUCTIVITIES% %DIPOLES% %DSM%
om_assemble -H2EM %GEOMETRY% %CONDUCTIVITIES% %EEG_ELECTRODES% %H2EM%
om_gain -EEG %HMINV% %DSM% %H2EM% %EEG_LEADFIELD%

rem Compute MEG gain matrix
om_assemble -H2MM %GEOMETRY% %CONDUCTIVITIES% %SQUIDS% %H2MM%
om_assemble -DS2MM %DIPOLES% %SQUIDS% %DS2MEG%
om_gain -MEG %HMINV% %DSM% %H2MM% %DS2MEG% %MEG_LEADFIELD%

rem Compute EIT gain matrix
om_assemble -EITSM %GEOMETRY% %CONDUCTIVITIES% %EEG_ELECTRODES% %EITSM%
om_gain -EEG %HMINV% %EITSM% %H2EM% %EIT_LEADFIELD%

rem Compute Internal Potential
om_assemble -H2IPM %GEOMETRY% %CONDUCTIVITIES% %INTERNAL_ELECTRODES% %IPHM%
om_assemble -DS2IPM %GEOMETRY% %CONDUCTIVITIES% %DIPOLES% %INTERNAL_ELECTRODES% %IPSM%
om_gain -IP %HMINV% %DSM% %IPHM% %IPSM% %IP_LEADFIELD%

