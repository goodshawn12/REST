function nut_loadVOI_Callback(VOIname,alsoview)
% NUT_LOADVOI_CALLBACK
%
% Load pre-defined VOI.
%
% Uses global NUTS and st
% VOIname is full path name of voxels.mat file 
% alsoview is option to view the loaded VOI or not (1=yes, 0=no)

global nuts st;
if ~exist('alsoview')
    alsoview=1;
end


if exist('VOIname')   %% for batching
    load(VOIname,'VOIvoxels');
else  %% for interactive gui use
	[VOIfilename, VOIpathname]=uigetfile('*.mat','Select VOI...(press cancel to use voxels already in current nuts)','voxels.mat');
	if ~(isequal(VOIfilename,0)|isequal(VOIpathname,0))
        load(fullfile(VOIpathname,VOIfilename),'VOIvoxels');
        nuts.VOIvoxels = VOIvoxels;
    else
        VOIvoxels=nuts.VOIvoxels;
	end
end


if alsoview
    nut_refresh_image;
    
    % convert from MEG coords to MRI coords, and make positive by
    % substracting st.bb
    VOIvoxels = nut_meg2mri(double(VOIvoxels));  %not needed now??
    %VOIvoxels = nut_mri2meg(VOIvoxels);  %not needed now??
    blobcenter = mean(VOIvoxels);
    translation_tfm = [ 1 0 0 -abs(st.bb(1,1))
                        0 1 0 -abs(st.bb(1,2))
                        0 0 1 -abs(st.bb(1,3))
                        0 0 0                1 ];
    VOIvoxels = nut_coordtfm(VOIvoxels,inv(translation_tfm));

    spm_image('init',st.vols{1}.fname); % reload structural MRI; clears VOI emphasis from nut_view_scan_region.m
    nut_spmfig_setup;

  	spm_orthviews('rmblobs',1)
    if exist('VOIvoxels','var')
        spm_orthviews('Reposition',blobcenter);

        %  so that the coloured blobs are not opaque, a contrast within the blob
		%  needs to exist; hence arbitrarily setting the last voxel to value=4
		%  while the rest remain value=1
        spm_orthviews('addcolouredblobs',1,VOIvoxels',ones(size(VOIvoxels,1),1),translation_tfm,[.3 1 .2]);
        st.vols{1}.blobs{1}.max=4;  % force desired color scale
        st.vols{1}.blobs{1}.min=0;
	    spm_orthviews('redraw');
    else
        warning('VOIvoxels not within this voxels.mat file or not created yet!');
    end
end

%nuts.VOIvoxels=VOIvoxels;  
clear VOIvoxels

% remove any preloaded lead field since we've obviously made it invalid by loading a new VOI
if(isfield(nuts,'Lp'))
    nuts = rmfield(nuts,'Lp');
end

nut_enabler;