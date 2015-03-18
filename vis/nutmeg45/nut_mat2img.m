function nut_mat2img(timept,outfile)
% nut_mat2img   Function to export activations to the Analyze (.img) format.
%               nut_timef_viewer must be open and displaying dataset when called.
%
% nut_mat2img(timept,outIMGfile)  [both inputs are optional]
%


global defaults rivets st beam
if ( ~isfield(rivets,'fig') || ~ishandle(rivets.fig) )
    help nut_mat2img
    return
end
if rivets.displaymode==2
    errordlg('Cannot export in 3D rendered mode. Please uncheck the "Display 3D" box first.')
    return
end

if(strcmp(spm('ver'),'SPM2'))
    if defaults.analyze.flip
        warndlg('WARNING: Your SPM default settings assume that MRIs have left and right sides flipped. This setting is incompatible with NUTMEG! Please type "edit spm_defaults" in the MATLAB command window, and change the line "defaults.analyze.flip=1" to "defaults.analyze.flip=0".')
    end
    defaults.analyze.flip=0; % needed to placate SPM2's trigger-happy flipping
elseif(strncmp(spm('ver'),'SPM8',4))
    if spm_flip_analyze_images
        errordlg('Your SPM default settings assume that MRIs have left and right sides flipped. This setting is incompatible with NUTMEG! Please type "edit spm_flip_analyze_images.m" in the MATLAB command window, and change the line "flip=1" to "flip=0".')
    end
else
    warning('not sure what to do about defaults.analyze.flip for this version of SPM');
end


% Get *.img filename
%-------------------
if nargin<2
    [outfile,outpath] = uiputfile([rivets.beamfilename '.img'],'Export as Analyze...');
    outfile=[outpath outfile];
    if ~ischar(outfile), return, end
end

% Prepare functional data
%------------------------
if rivets.timefmode
    if size(rivets.s,2)==1
        timeopt='selected';
        imgvec=rivets.s(:,1,rivets.freqselect);
        select = any(rivets.threshold(:,1,:),3);
    elseif nargin>0
        timeopt='selected';
        imgvec=rivets.s(:,timept,rivets.freqselect);
        select = any(rivets.threshold(:,timept,:),3);
    else
        timeopt=questdlg('Which time window(s) would you like to export?','','Selected','Area under curve of all','Selected');
        switch timeopt
            case 'Selected'
                imgvec=rivets.s(:,rivets.timeselect,rivets.freqselect);
                select = any(rivets.threshold(:,rivets.timeselect,:),3);
            case 'AUC of all'
                imgvec=trapz(beam.timepts,rivets.s(:,:,rivets.freqselect),2);
                select = any(rivets.threshold(:,rivets.timeselect,:),3);
            otherwise
                return
        end
    end
else
    if nargin>0
        timeopt='selected';
        imgvec=rivets.s(:,timept);
        select = any(rivets.threshold(:,timept,:),3);
    else
        handles=guihandles(rivets.fig);
        time = str2double(get(handles.nut_time_text,'String'));
        time2 = str2double(get(handles.nut_time2_text,'String'));
        timept = dsearchn(beam.timepts,time);   % find nearest index in timevector to given time
        timept2 = dsearchn(beam.timepts,time2);
        if (timept2-timept==0)
            timeopt='selected';
            imgvec=rivets.s(:,timept);
            select = any(rivets.threshold(:,timept,:),3);
        else
            timeopt='area under curve of selected';
            imgvec=trapz(beam.timepts(timept:timept2),rivets.s(:,timept:timept2),2);
            select = any(rivets.threshold(:,rivets.timeselect,:),3);
        end
    end
end

%if ~any(select), errordlg('No activations were above threshold.'), return, end
imgvec(~select)=0;

% Co-align with structural MRI
%-----------------------------
data = nut_vector2vol(imgvec,beam);

Vm.dim=size(data);
Vm.fname=outfile;
Vm.pinfo=[1;0;0];
Vm.mat=st.vols{1}.blobs{1}.mat;
% datatype handled differently between SPM2 and later versions
% 64 means double precision
if(strcmp(spm('ver'),'SPM2'))
    Vm.dim(4) = 64;
else
    Vm.dt = [64 0];
end

% fuck spm -- if true, removes NaN values for better compatibility with external image viewers
fuckspm = true;
if(fuckspm)
    data(isnan(data))=0;   % for compatibility with external programs
end

% write output
spm_write_vol(Vm,data);

% reslice if necessary
answer=questdlg('Some programs may not display the overlay correctly over the structural MRI. If this is the case, try to reslice. Do you want to do this now?','NUTMEG question','Yes','No','No');
if strcmp(answer,'Yes')
    V=spm_vol(beam.coreg.mripath);
    V(2)=spm_vol(outfile);
    answer=questdlg('Do you want to reslice other structural MRI files to the same coordinates (e.g., skull stripped MRI)?','NUTMEG question','Yes','No','Yes');
    if strcmp(answer,'Yes')
        [fo, po] = uigetfile('*.img;*.nii','Select other MRIs to reslice','MultiSelect','on');
        if ~iscell(fo), fo={fo}; end
        for k=1:length(fo)
            V(2+k)=spm_vol(fullfile(po,fo{k}));
        end
    end
    flags = struct('which',2,'mean',0);
    spm_reslice(V,flags)
    [outpath,outfile,ext]=fileparts(outfile);
    outfile = fullfile(outpath,['r' outfile],ext);
end 

if nargin==0
    msgbox(sprintf('The functional data of the %s time window(s) was exported to Analyze format as %s.',lower(timeopt),outfile))
end
    
%if ( mrifile(1)=='w' || strcmp(mrifile,'T1') || strncmp(mrifile,'avg',3) )
% if blobs are normalized, we have a quick and easy job
    
% Not sure this ever worked correctly, it does not seem to work now...
%     Vs=spm_vol(beam.coreg.mripath);
% 
%     spixdim=nut_pixdim(Vs);
% 
%     Vm.fname=outfile;
%     xdim=floor(Vs.dim(1)*abs(spixdim(1))/beam.voxelsize(1));
%     ydim=floor(Vs.dim(2)*abs(spixdim(2))/beam.voxelsize(2));
%     zdim=floor(Vs.dim(3)*abs(spixdim(3))/beam.voxelsize(3));
%     Vm.dim = [xdim ydim zdim];
% 
%     % datatype handled differently between SPM2 and later versions
%     % 64 means double precision
%     if(strcmp(spm('ver'),'SPM2'))
%         Vm.dim(4) = 64;
%     else
%         Vm.dt = [64 0];
%     end
%  
%     Vs_origmat=[spixdim(1) 0 0 -(Vs.dim(1)+1)/2*spixdim(1); 0 spixdim(2) 0 -(Vs.dim(2)+1)/2*spixdim(2); 0 0 spixdim(3) -(Vs.dim(3)+1)/2*spixdim(3); 0 0 0 1];
%     Vs_rot=Vs.mat*inv(Vs_origmat);
% 
% %     Vm.mat=[beam.voxelsize(1) 0 0 -(xdim+1)/2*beam.voxelsize(1);0 beam.voxelsize(2) 0 -(ydim+1)/2*beam.voxelsize(2);0 0 beam.voxelsize(3) -(zdim+1)/2*beam.voxelsize(3);0 0 0 1];
%     Vm.mat=[beam.voxelsize(1) 0 0 -xdim/2*beam.voxelsize(1)-30;0 beam.voxelsize(2) 0 -ydim/2*beam.voxelsize(2)-30;0 0 beam.voxelsize(3) -zdim/2*beam.voxelsize(3)-30;0 0 0 1];
%     Vm.mat=Vs_rot*Vm.mat;
%     Vm.pinfo=Vs.pinfo;
%     Vm.descrip=[];
% 
%     iVmmat=inv(Vm.mat);
%     mricoord=nut_coordtfm(rivets.voxelsMRI,iVmmat);
% 
%     % Avoid errrors due to functional blobs that are outside structural MRI
%     insidebrain=find( all(mricoord>0,2) & mricoord(:,1)<=xdim & mricoord(:,2)<=ydim & mricoord(:,3)<=zdim ); % use only voxels that are inside structural MRI
%     mricoord=mricoord(insidebrain,:); 
%     imgvec=imgvec(insidebrain);
%     clear insidebrain

%     global beam st
%     Vm.fname=outfile;
%     Vm.mat=st.vols{1}.blobs{1}.mat;
%     Vm.dim=size(st.vols{1}.blobs{1}.vol);
%     Vm.mat(1:3,4)=Vm.mat(1:3,4)-7;
%     % datatype handled differently between SPM2 and later versions
%     % 16 means float, 64 means double precision
%     if(strcmp(spm('ver'),'SPM2'))
%         Vm.dim(4) = 64;
%     else
%         Vm.dt = [64 0];
%     end
%     Vm.pinfo=[1;0;0];
% 
%     data = st.vols{1}.blobs{1}.vol;

%for jj=1:zdim,
%    Vspm=spm_write_plane(Vm,squeeze(data(:,:,jj)),jj);
%end
%spm_close_vol(Vout);


