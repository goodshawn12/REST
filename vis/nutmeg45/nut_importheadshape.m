function nut_importheadshape(headshapefile)
% Read the Head Shape file, get the headshape points and the fiducials       

global  coreg st screws

if(~exist('headshapefile','var'))
    [headshapefile, headshapepath]=uigetfile({'*.ahs;*.hsp;*.pos;*.shape;*.mat;*.elc;hs_file.*'},'Please select the head shape.');
    if isequal(headshapepath,0)|isequal(headshapefile,0)
        % screw this, user hit cancel
        return;
    else
        headshapefile = fullfile(headshapepath,headshapefile);
    end
end

[headshapepath,hsfile_no_ext,hstype] = fileparts(headshapefile);

% BTi/4D headshape file "hs_file" -- requires Eugene Kronberg's PDF4D toolbox (biomag.wikidot.com)
if(strcmp(hsfile_no_ext,'hs_file'))
    tmp = read_hs_file(headshapefile);
    fiducials_meg=[tmp.index.lpa tmp.index.rpa tmp.index.nasion tmp.index.cz tmp.index.inion]'*1000; % convert from m to mm
    % also present:
    % tmp.index.cz;
    % tmp.index.inion;
    coreg.hsCoord=tmp.point'*1000 % convert from m to mm
    coreg.hsCoord=[coreg.hsCoord;fiducials_meg];
else
    switch(hstype)
        case '.hsp'
            [crap1,fid_x,fid_y,fid_z] = textread(headshapefile,'%s%f%f%f',3,'headerlines',9);
            fiducials_meg=[fid_x,fid_y,fid_z]*1000;  % convert from m to mm

            [hs_x,hs_y,hs_z] = textread(headshapefile,'%f%f%f','headerlines',14);
            coreg.hsCoord=[hs_x,hs_y,hs_z]*1000; % convert from m to mm

            clear crap1 crap2;
        case '.ahs'
            [a,b,c,d,e,f,g,h,i,j,k,l,m]=textread(headshapefile,'%f%f%f%n%n%n%n%n%n%n%n%s%s',2, 'headerlines',1);
            [a(3),b(3),c(3),d,e,f,g,h,i,j,k,l]=textread(headshapefile,'%f%f%f%n%n%n%n%n%n%n%n%s',1, 'headerlines',3);
            fiducials_meg=[a,b,c]*10;  % convert from cm to mm
            fiducials_meg=fiducials_meg';
            [a1,b1,c1,d,e,f,g,h,i,j,k,l]=textread(headshapefile,'%f%f%f%n%n%n%n%n%n%n%n%s','headerlines',6);
            coreg.hsCoord=[a1,b1,c1]*10;  % convert from cm to mm
            coreg.hsCoord=[nuts.coreg.hsCoord;fiducials_meg'];
            clear a1 b1 c1 a b c d e f g h i j k l m;
        case '.pos'
            [num] = textread(headshapefile,'%f',1);
            [ind,hs_x,hs_y,hs_z] = textread(headshapefile,'%s%f%f%f',num,'headerlines',1);
            [label,fid_x,fid_y,fid_z] = textread(headshapefile,'%s%f%f%f',3,'headerlines',1+num);
            fiducials_meg=[fid_x,fid_y,fid_z]'*10; % convert from cm to mm
            coreg.hsCoord=[hs_x,hs_y,hs_z]*10; % convert from cm to mm
            coreg.hsCoord=[coreg.hsCoord;fiducials_meg'];
        case '.shape'
            [num] = textread(headshapefile,'%f',1);
            [hs_x,hs_y,hs_z] = textread(headshapefile,'%f%f%f',num,'headerlines',1);
            coreg.hsCoord=[hs_x,hs_y,hs_z]*10; % convert from cm to mm
        case '.elc' % Polaris digitizer electrode positions
            ftelec = ft_read_sens(headshapefile);
            coreg.hsCoord = ftelec.elecpos;
        otherwise
            try  % could be EEG electrode locations; use FieldTrip to import
                ftelec = ft_read_sens(headshapefile);
                coreg.hsCoord = ftelec.elecpos;
            catch
                errordlg('Unsupported headshape format.');
                return
            end
    end
end


voxelsize = [4 4 4];
shift = [ (min(coreg.hsCoord(:,1)) + max(coreg.hsCoord(:,1)) - 256)*0.5
          (min(coreg.hsCoord(:,2)) + max(coreg.hsCoord(:,2)) - 256)*0.5
          (min(coreg.hsCoord(:,3)) + max(coreg.hsCoord(:,3)) - 256)*0.5 ];

% translation_tfm shifts points to be in center of created MRI volume and sets voxelsize
translation_tfm = [ voxelsize(1)            0            0 shift(1,:)
                               0 voxelsize(2)            0 shift(2,:)
                               0            0 voxelsize(3) shift(3,:)
                               0            0            0          1 ];

if(strfind(coreg.mripath,'nutmeg/templates/blank')) % if MRI loaded, superimpose headshape points
    hsvoxels = nut_coordtfm(coreg.hsCoord,inv(translation_tfm));

    imgformat = 2; % Analyze format, type 2 (i.e., uint8);
	maxvalue = spm_type(imgformat,'maxval');
	hsimg = zeros(256/voxelsize(1),256/voxelsize(2),256/voxelsize(3));
	for i=1:length(hsvoxels)
        hsimg(ceil(hsvoxels(i,1)),ceil(hsvoxels(i,2)),ceil(hsvoxels(i,3)))=maxvalue;
        hsimg(ceil(hsvoxels(i,1)),ceil(hsvoxels(i,2)),floor(hsvoxels(i,3)))=maxvalue;
        hsimg(ceil(hsvoxels(i,1)),floor(hsvoxels(i,2)),ceil(hsvoxels(i,3)))=maxvalue;
        hsimg(floor(hsvoxels(i,1)),ceil(hsvoxels(i,2)),ceil(hsvoxels(i,3)))=maxvalue;
        hsimg(floor(hsvoxels(i,1)),floor(hsvoxels(i,2)),ceil(hsvoxels(i,3)))=maxvalue;
        hsimg(floor(hsvoxels(i,1)),ceil(hsvoxels(i,2)),floor(hsvoxels(i,3)))=maxvalue;
        hsimg(ceil(hsvoxels(i,1)),floor(hsvoxels(i,2)),floor(hsvoxels(i,3)))=maxvalue;
        hsimg(floor(hsvoxels(i,1)),floor(hsvoxels(i,2)),floor(hsvoxels(i,3)))=maxvalue;
	end
	
	hsvol.fname = fullfile(headshapepath,[hsfile_no_ext 'shape.img']);
	hsvol.dim = [size(hsimg,1) size(hsimg,2) size(hsimg,3)];

    if(strcmp(spm('ver'),'SPM2'))
        Vm.dim(4) = imgformat;
    else
        Vm.dt = [imgformat 0];
    end

    hsvol.mat = translation_tfm;
	hsvol.pinfo = [1;0;0];
	hsvol=spm_write_vol(hsvol,hsimg);
	
    disp('If I only had a brain...');
    coreg.mripath = hsvol.fname;
    nut_load_image(coreg.mripath);
    
    coreg.fiducials_mri_mm = fiducials_meg(1:3,:);
    
    set(screws.coreg.handles.nut_nose_text,'String',sprintf('MRI: %.1f %.1f %.1f',coreg.fiducials_mri_mm(1,:)));
    set(screws.coreg.handles.nut_left_text,'String',sprintf('MRI: %.1f %.1f %.1f',coreg.fiducials_mri_mm(2,:)));
    set(screws.coreg.handles.nut_right_text,'String',sprintf('MRI: %.1f %.1f %.1f',coreg.fiducials_mri_mm(3,:)));
else  % if no MRI loaded, create Analyze format image volume containing headshape points
    hs_blobs = nut_coordtfm(coreg.hsCoord,inv(translation_tfm));
    spm_orthviews('rmblobs',1);
	spm_orthviews('addblobs',1,hs_blobs',zeros(size(hs_blobs,1),1),coreg.meg2mri_tfm*translation_tfm);
    if(strcmp(spm('ver'),'SPM2'))
        % only for SPM2; SPM8 freaks out if we delete the colorbar
        delete(st.vols{1}.blobs{1}.cbar);
    end
	spm_orthviews('redraw');
end

nut_coreg_enabler



