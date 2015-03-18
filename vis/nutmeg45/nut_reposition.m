function nut_reposition(pos)
% calls spm_orthviews('Reposition',pos), then resyncs MNI and MEG coords

spm_orthviews('Reposition',pos);
nut_image('shopos');