function nut_mni_fiducials_4spm8(outputdir)

global nuts
if isfield(nuts.coreg,'norm_mripath')
    nuts.coreg.fiducials_mni_mm=nut_mri2mni(nuts.coreg.fiducials_mri_mm,nuts.coreg,1);
else
    error('you must first specificy a normalized MRI');
end

lpa=nuts.coreg.fiducials_mni_mm(1,:);
rpa=nuts.coreg.fiducials_mni_mm(2,:);
nas=nuts.coreg.fiducials_mni_mm(3,:);

save([outputdir 'fidmni_lpa.mat'],'lpa')
save([outputdir 'fidmni_rpa.mat'],'rpa')
save([outputdir 'fidmni_nas.mat'],'nas')
