function nut_ftmriN2coreg(mripath,mristruct,fid,mmvoxflag)
% nut_ftmriN2coreg(mripath,mristruct,fid,mmvoxflag)
%
% populates the global 'coreg' structure in Nutmeg with the MRI and
% fiducial information
%
% mripath:   EITHER existing file fullpath name, both .nii or .img supported
%            OR desired fullpath name Plus 'mristruct' 
%            OR [] empty, default T1 MNI brain set
%              If you start with CTF .mri file, see nutmeg wiki for steps to
%              correctly convert to Analyze
% mristruct: FT structure from ft_read_mri (enter [] if mripath file already exists)
%
% fid: generically can be .txt file with list of fiducials 
%      IF in mm, then set mmvoxflag=1
%      ELSEIF in voxels, then set mmvoxflag=2
%      OR can be *.hdm file from CTF (set mmvoxflag=1)
%      see nut_load_fiducials.m for more details
%
% mmvoxflag = 1 if fiducials are in millimeters, =2 if in voxels
%


global coreg

if exist(mripath,'file')
    if (strcmp(ft_filetype(mripath),'nifti') || strcmp(ft_filetype(mripath),'analyze_img'))
        coreg.mripath=mripath;
    else
        error('please convert your MRI file to Analyze or Nifti first, ensuring fiducial info still correct');
    end        
    defmri=0;
elseif ~isempty(mripath)
    if isstruct(mristruct)
        cfg1 = [];
        cfg1.parameter   = 'anatomy';
        cfg1.filename    = mripath;
        cfg1.filetype    = 'nifti'; % or 'nifti_img' if you still want img/hdr pair
        cfg1.coordsys    = 'spm'; % 'ctf' also works for Yokogawa...does it generically work?
        ft_volumewrite(cfg1, mristruct)
%         mri=ft_read_mri(mripath);  % or 'mrinifti.img'
        coreg.mripath=mripath;
    else
        error('you must include mristruct if mripath not exist')
    end
    defmri=0;
else
    coreg.mripath = which('wNormT1neuroax.img');
    coreg.brainrender_path = which('render_wNormT1neuroax_seg1.mat');
    warning('Default T1 MNI brain used. Set fiducial coregistration correctly in GUI')
    % FIXME what about fiducials in this case?  
    defmri=1;
end
if isempty(fileparts(coreg.mripath))
    error('best to specify full path name of MRI');
end

if defmri==0
    V=spm_vol(coreg.mripath);
    if mmvoxflag==1
        nut_load_fiducials(fid,V); % coreg is already assumed to be global here
    elseif mmvoxflag==2
        fidvox=load(fid);
        fiducials_mri_mm=nut_coordtfm(fidvox,V.mat);
        save('fiducialsmm.txt','-ascii','fiducials_mri_mm');
        nut_load_fiducials('fiducialsmm.txt',V); % coreg is already assumed to be global here
    else
        error('mmvoxflag must be 1 or 2')
    end
else
    coreg.meg2mri_tfm=eye(4);
end
    
mri=ft_read_mri(coreg.mripath);
% interestingly, this line below used to work but does not anymore.
if isfield(mri,'transform')
    coreg.v2_meg2mri_tfm = mri.transform/mri.hdr.mat;
end;

