%%
% Projects the EEG sensor coordinates to the MR headshape for cases where
% the MR is available for the subject. Otherwise, allows warping of a
% template MR to the EEG sensor and fiducial coordinates (only for EEG
% systems with fiducial coregistration).
% @author Daniel D.E. Wong
%%
function nut_eegmralign(fitType)

global nuts st

if ~exist('PrepareTriangleMesh.m','file')
    errordlg('This function requires the Helsinki BEM library in your Matlab path. It can be downloaded from http://peili.hut.fi/BEM')
    return
end

if nargin<1 
    if isfield(nuts.meg,'lpa')
        fitType = questdlg('Select operation:','Coregistration Operations','Project Sensors','Warp MR','Cancel','Project Sensors');
    else    % For EEG systems without fiducial registration, the Warp option is not available.
        fitType = 'Project Sensors';
    end
end
if strcmp(fitType,'Cancel'); return; end;

% Load MRI data
if strcmp(fitType,'Project Sensors')
    Y = spm_read_vols(st.vols{1});
else
    [Y,XYZ] = spm_read_vols(st.vols{1});
    XYZ = XYZ';
end

% Create mrifid data structure
if ~isfield(nuts.coreg,'volfile')
    [file,path] = uigetfile('*_vol.mat','Open Mesh File');
    if isequal(file,0); errordlg('No mesh file loaded'); return; end;
    nuts.coreg.volfile = fullfile(path,file);
end
load(nuts.coreg.volfile);
% Convert to MEG coordinates
for i = 1:length(vol.bnd)
    % Handle fieldtrip vols
    if isfield(vol.bnd(i),'tri'); vol.bnd(i).faces = vol.bnd(i).tri; end;
    if isfield(vol.bnd(i),'pnt'); vol.bnd(i).vertices = vol.bnd(i).pnt; end;
    
    vol.bnd(i).vertices = nut_mri2meg(vol.bnd(i).vertices);
end
while 1
    answer = questdlg(['Increase scalp density? (' num2str(size(vol.bnd(1).faces,1),1) 'faces)'],'Display Mesh','Yes','No','No');
    if strcmp(answer,'Yes')
        vol.bnd(1) = mesh_refine_tri4(vol.bnd(1));
    else
        break;
    end
end
mrifid.hsp = PrepareTriangleMesh(vol.bnd(1).vertices,vol.bnd(1).faces(:,[1 3 2]));
mrifid.lpa = nut_mri2meg(nuts.coreg.fiducials_mri_mm(1,:));
mrifid.rpa = nut_mri2meg(nuts.coreg.fiducials_mri_mm(2,:));
mrifid.nasion = nut_mri2meg(nuts.coreg.fiducials_mri_mm(3,:));

% Create eegfid data structure
eegfid.sensorCoord = nuts.meg.sensorCoord;
if strcmp(fitType,'Warp MR') % nuts.meg.lpa is not always available and does not seem to be necessary for 'project sensors'
    eegfid.lpa = nuts.meg.lpa;
    eegfid.rpa = nuts.meg.rpa;
    eegfid.nasion = nuts.meg.nasion;
    if isfield(nuts.meg,'inion'); eegfid.inion = nuts.meg.inion; mrifid.inion = [-90 0 25]; end;    % MRI inion is for standard_mri.img
end

if strcmp(fitType,'Project Sensors') | strcmp(fitType,'Warp MR')
    disp('Projecting sensors to MR head surface...');
    
    if ~exist('vol','var')
        mrSensorCoord = zeros(size(eegfid.sensorCoord));
        for i = 1:size(eegfid.sensorCoord,1)
            mrSensorCoord(i,:) = mrifid.hsp(dsearchn(mrifid.hsp,eegfid.sensorCoord(i,:)),:);
        end
    else
        mrSensorCoord = nuts.meg.sensorCoord;
        % Some EEG caps go further down that the scalp surface, so we
        % cannot project these sensors
        for i = find( nuts.meg.sensorCoord(:,3) >= min(mrifid.hsp.mp(:,3)) )' %1:size(nuts.meg.sensorCoord,1)
            senstrimap = dsearchn(mrifid.hsp.mp,nuts.meg.sensorCoord(i,:));
            mrSensorCoord(i,:) = nuts.meg.sensorCoord(i,:) - mrifid.hsp.n(senstrimap,:).*dot(mrifid.hsp.n(senstrimap,:),nuts.meg.sensorCoord(i,:)-mrifid.hsp.mp(senstrimap,:))/sum(mrifid.hsp.n(senstrimap,:).^2);
        end
    end
    
    % Save sensor locations
    if strcmp(fitType,'Project Sensors')
        nuts.meg.sensorCoord = mrSensorCoord;
        nuts.meg.lpa = nut_mri2meg(nuts.coreg.fiducials_mri_mm(1,:));
        nuts.meg.rpa = nut_mri2meg(nuts.coreg.fiducials_mri_mm(2,:));
        nuts.meg.nasion = nut_mri2meg(nuts.coreg.fiducials_mri_mm(3,:));
    end
end
if strcmp(fitType,'Warp MR')
    % Warp template MR headshape by mapping projected sensors to actual
    % sensors (this provides an approximate interpolated headshape)
    disp('Interpolating headshape...');
    
    mrSensorCoord = mrSensorCoord(nuts.meg.goodchannels,:);             % Only goodchannels should be used for head surface interpolation
    eegfid.sensorCoord = eegfid.sensorCoord(nuts.meg.goodchannels,:);
    warp_vertices = tps_warp([mrifid.nasion; mrifid.lpa; mrifid.rpa; mrSensorCoord],[eegfid.nasion; eegfid.lpa; eegfid.rpa; eegfid.sensorCoord],mrifid.hsp.p);
    %warp_vertices = tps_warp(mrSensorCoord,eegfid.sensorCoord,mrifid.hsp.p);
    eegfid.hsp = PrepareTriangleMesh(warp_vertices,vol.bnd(1).faces(:,[1 3 2]));
    %figure, trimesh(eegfid.hsp.e,eegfid.hsp.p(:,1),eegfid.hsp.p(:,2),eegfid.hsp.p(:,3)); hold on; plot3(nuts.meg.sensorCoord(nuts.meg.goodchannels,1),nuts.meg.sensorCoord(nuts.meg.goodchannels,2),nuts.meg.sensorCoord(nuts.meg.goodchannels,3),'o'); hold off; return;
    
    % Obtain 10-20 sensor locations from interpolated headshape
    disp('Generating 10-10 sensor locations...');
    refpointsEEG = find1020Refs(eegfid); %refpointsEEG = refpointsEEG(5:end,:);    % Exclude fiducials
    %figure, trimesh(eegfid.hsp.e,eegfid.hsp.p(:,1),eegfid.hsp.p(:,2),eegfid.hsp.p(:,3)); hold on; plot3(refpointsEEG(:,1),refpointsEEG(:,2),refpointsEEG(:,3),'o'); hold off; return;
    
    % Obtain 10-20 sensor locations from original MR headshape
    refpointsMR = find1020Refs(mrifid); %refpointsMR = refpointsMR(5:end,:);
    %figure, trimesh(mrifid.hsp.e,mrifid.hsp.p(:,1),mrifid.hsp.p(:,2),mrifid.hsp.p(:,3)); hold on; plot3(refpointsMR(:,1),refpointsMR(:,2),refpointsMR(:,3),'o'); hold off; return;
    
    % Warp original MR by mapping 10-20 sensor locations
    disp('Computing interpolation splines...')
    d = [3 3 3 0 0 0];
    c = spm_bsplinc(st.vols{1},d);
    
    disp('Warping MR...')
    XYZ = single(nut_mri2meg(XYZ));
    XYZ = single(tps_warp(single(refpointsEEG), single(refpointsMR), XYZ));
    XYZ = double(nut_mm2voxels(nut_meg2mri(XYZ)));
    
    disp('Interpolating...')
    Y = spm_bsplins(c,XYZ(:,1),XYZ(:,2),XYZ(:,3),d);
    Y = reshape(Y,st.vols{1}.dim(1:3));
    size(Y)
    clear XYZ;
    
    % Save warped MR
    disp('Writing warped MR...')
    V = st.vols{1};
    [file,path] = uiputfile('*.img','Save Warped MRI');
    if ~isequal(file,0) & ~isequal(path,0)
        V.fname=strcat(path,file);
        spm_write_vol(V,Y);
    end
    
    % Warp BEM points
    [file,path] = uigetfile('*_vol.mat','Open Source Mesh File for Warping');
    if ~isequal(file,0) & ~isequal(path,0)
        load(strcat(path,file));
        for i = 1:length(vol.bnd)
            vol.bnd(i).vertices = nut_mri2meg(vol.bnd(i).vertices);
            vol.bnd(i).vertices = tps_warp(refpointsMR, refpointsEEG, vol.bnd(i).vertices);
            vol.bnd(i).vertices = nut_meg2mri(vol.bnd(i).vertices);
            %vol.bnd(i).vertices = vol.bnd(i).vertices-repmat(minXYZ,size(vol.bnd(i).vertices,1),1);
            %vol.bnd(i).vertices = vol.bnd(i).vertices-repmat((dim'+1)/2,size(vol.bnd(i).vertices,1),1);
        end
        [file,path] = uiputfile('*_vol.mat','Save Warped Mesh');  % BEM points are saved in MRI coordinates
        if ~isequal(file,0) & ~isequal(path,0)
            save(strcat(path,file),'vol');
        end
    end
    
    % Warp fiducial coordinates
    fids = nut_mri2meg(nuts.coreg.fiducials_mri_mm);
    fids = tps_warp(refpointsMR, refpointsEEG, fids);
    fids = double(nut_meg2mri(fids));
    %fids = fids-repmat(minXYZ,size(fids,1),1);      % Re-adjust fids to correspond to new MRI coordinate origin
    %fids = double(fids-repmat((dim'+1)/2,size(fids,1),1));
    [file,path] = uiputfile('*.txt','Save Warped Fiducial Points');
    if ~isequal(file,0) & ~isequal(path,0)
        save(strcat(path,file),'fids','-ascii');
    end
end



% ************************
% * Supporting Functions *
% ************************



%%
% Finds 10-20 sensor coordinates.
% @param XYZ the MR headshape points.
% @param fid the fiducial structure.
% @return refpoints reference points in the 10-20 system as a 25x3 matrix. Refer
% to comments to determine which index corresponds to what location.
%%
function refpoints = find1020Refs(fid)

r = 150;                                                                % Search radius

refpoints = zeros(101,3);
refpoints(1,:) = fid.lpa;                                               % lpa
refpoints(2,:) = fid.rpa;                                               % rpa
refpoints(3,:) = fid.nasion;                                            % nasion

% Find inion
if ~isfield(fid,'inion')                                                % inion
    tri_idx = dsearchn(fid.hsp.mp,[-r 0 0]);
    refpoints(4,:) = triangle_intersect(fid.hsp.p(fid.hsp.e(tri_idx,:),:), [-r 0 0]);   % autoselect inion
else
    refpoints(4,:) = fid.inion;                                         
end
refpoints(4,2) = 0;                                                     % Ensure inion is in xz-plane

% Obtain saggital reference points
radpts = zeros(181,3);      % 0-180 deg
radpts(1,:) = fid.nasion;
deg = [0:180]'; pt = [r*cosd(deg), zeros(181,1), r*sind(deg)];
tri_idx = dsearchn(fid.hsp.mp,pt);
for deg = 0:180
    radpts(deg+1,:) = triangle_intersect(fid.hsp.p(fid.hsp.e(tri_idx(deg+1),:),:), pt(deg+1,:));
end
%figure, plot3(radpts(:,1),radpts(:,2),radpts(:,3),'.'); hold on; trimesh(fid.hsp.e,fid.hsp.p(:,1),fid.hsp.p(:,2),fid.hsp.p(:,3)); hold off; return

inionidx = dsearchn(radpts,refpoints(4,:));
circ = cumsum(sqrt(sum(diff(radpts).^2,2)));
for deg = 0:179
    fracpos = circ(deg+1)/circ(inionidx-1);
    if fracpos <= 0.1; refpoints(5,:) = radpts(deg+1,:);             % FPz
    elseif fracpos <= 0.2; refpoints(6,:) = radpts(deg+1,:);         % AFz
    elseif fracpos <= 0.3; refpoints(7,:) = radpts(deg+1,:);         % Fz
    elseif fracpos <= 0.4; refpoints(8,:) = radpts(deg+1,:);         % FCz
    elseif fracpos <= 0.5; refpoints(9,:) = radpts(deg+1,:);         % Cz
    elseif fracpos <= 0.6; refpoints(10,:) = radpts(deg+1,:);        % CPz
    elseif fracpos <= 0.7; refpoints(11,:) = radpts(deg+1,:);        % Pz
    elseif fracpos <= 0.8; refpoints(12,:) = radpts(deg+1,:);        % POz
    elseif fracpos <= 0.9; refpoints(13,:) = radpts(deg+1,:);        % Oz
    else; break; end;
end
%figure, plot3(refpoints(5:13,1),refpoints(5:13,2),refpoints(5:13,3),'o'); hold on; trimesh(fid.hsp.e,fid.hsp.p(:,1),fid.hsp.p(:,2),fid.hsp.p(:,3)); hold off; return

% Obtain circumferential reference points
ori = (refpoints(5,:) + refpoints(13,:))/2;     % Place origin between Fpz and Oz
planenorm = cross([0 1 0], refpoints(5,:)-ori); % Define circumferential plane
d = -dot(refpoints(5,:),planenorm);

radpts = zeros(361,3);
%radpts(1,:) = refpoints(5,:);   % First point is at Fpz
%radpts(181,:) = refpoints(13,:); % Mid point is at Oz
%radpts(361,:) = refpoints(5,:); % Last point is again at Fpz
deg = [0:360]'; pt = [r*cosd(deg), r*sind(deg), zeros(361,1)];
p_c = fid.hsp.p-repmat(ori,size(fid.hsp.p,1),1);
mp_c = fid.hsp.mp-repmat(ori,size(fid.hsp.mp,1),1);
tri_idx = dsearchn(mp_c,pt);
for deg = 0:360
    %if deg == 180 | deg == 360; continue; end;
    radpts(deg+1,:) = triangle_intersect(p_c(fid.hsp.e(tri_idx(deg+1),:),:),pt(deg+1,:)) + ori;
end
%figure, plot3(radpts(:,1),radpts(:,2),radpts(:,3),'.'); hold on; trimesh(fid.hsp.e,fid.hsp.p(:,1),fid.hsp.p(:,2),fid.hsp.p(:,3)); hold off; return

% Left hemisphere
circ = cumsum(sqrt(sum(diff(radpts).^2,2)));
for deg = 0:360
    fracpos = circ(deg+1)/circ(180);
    if fracpos <= 0.1; refpoints(14,:) = radpts(deg+1,:);              % Fp1
    elseif fracpos <= 0.2; refpoints(15,:) = radpts(deg+1,:);          % AF1
    elseif fracpos <= 0.3; refpoints(16,:) = radpts(deg+1,:);          % F7
    elseif fracpos <= 0.4; refpoints(17,:) = radpts(deg+1,:);          % FT7
    elseif fracpos <= 0.5; refpoints(18,:) = radpts(deg+1,:);          % T7
    elseif fracpos <= 0.6; refpoints(19,:) = radpts(deg+1,:);          % TP7
    elseif fracpos <= 0.7 refpoints(20,:) = radpts(deg+1,:);           % P7
    elseif fracpos <= 0.8; refpoints(21,:) = radpts(deg+1,:);          % PO7
    elseif fracpos <= 0.9; refpoints(22,:) = radpts(deg+1,:);          % O1
    else
        deg1 = deg;
        break;
    end
end

% Right hemisphere
for deg = deg1:360
    fracpos = (circ(deg+1)-circ(180))/(circ(end)-circ(180));
    if fracpos <= 0.1; refpoints(23,:) = radpts(deg+1,:);              % O2
    elseif fracpos <= 0.2; refpoints(24,:) = radpts(deg+1,:);          % PO8
    elseif fracpos <= 0.3; refpoints(25,:) = radpts(deg+1,:);          % P8
    elseif fracpos <= 0.4; refpoints(26,:) = radpts(deg+1,:);          % TP8
    elseif fracpos <= 0.5; refpoints(27,:) = radpts(deg+1,:);          % T8
    elseif fracpos <= 0.6; refpoints(28,:) = radpts(deg+1,:);          % FT8
    elseif fracpos <= 0.7; refpoints(29,:) = radpts(deg+1,:);          % F8
    elseif fracpos <= 0.8; refpoints(30,:) = radpts(deg+1,:);          % FP8
    elseif fracpos <= 0.9; refpoints(31,:) = radpts(deg+1,:);          % Fp2
    else; break; end;
end
%figure, plot3(refpoints(5:31,1),refpoints(5:31,2),refpoints(5:31,3),'o'); hold on; trimesh(fid.hsp.e,fid.hsp.p(:,1),fid.hsp.p(:,2),fid.hsp.p(:,3)); hold off; return

refpoints(32:41,:) = getEEG34coord(refpoints(15,:),refpoints(6,:),refpoints(30,:),r,fid);
refpoints(42:51,:) = getEEG34coord(refpoints(16,:),refpoints(7,:),refpoints(29,:),r,fid);
refpoints(52:61,:) = getEEG34coord(refpoints(17,:),refpoints(8,:),refpoints(28,:),r,fid);
refpoints(62:71,:) = getEEG34coord(refpoints(18,:),refpoints(9,:),refpoints(27,:),r,fid);
refpoints(72:81,:) = getEEG34coord(refpoints(19,:),refpoints(10,:),refpoints(26,:),r,fid);
refpoints(82:91,:) = getEEG34coord(refpoints(20,:),refpoints(11,:),refpoints(25,:),r,fid);
refpoints(92:101,:) = getEEG34coord(refpoints(21,:),refpoints(12,:),refpoints(24,:),r,fid);



%%
% Computes electrode coordinates of quarter-line electrodes (i.e. F3 and
% F4).
% @param XYZ headshape coordinates.
% @param x7 left circumferential coordinate (i.e. F7).
% @param xz center-line coordinate (i.e. Fz).
% @param x8 right circumferential coordinate (i.e. F8).
% @return x3 left quarter-line coordinate (i.e. F3).
% @return x4 right quarter-line coordinate (i.e. F4).
%%
function refpoints = getEEG34coord(x7,xz,x8,r,fid)

% Find points on a circle and project onto a plane
planevec1 = x8 - x7; planevec1 = planevec1/norm(planevec1);
planevec2 = xz - x7; planevec2 = planevec2/norm(planevec2);
planenorm = cross(planevec2,planevec1);     % Plane normal pointing forward
d = -dot(x7,planenorm);                     % Ensure the plane runs through X7

ori = (x7 + x8)/2;
p_a = fid.hsp.p-repmat(ori,size(fid.hsp.p,1),1);
mp_a = fid.hsp.mp-repmat(ori,size(fid.hsp.mp,1),1);
deg = [1:180]'; pt = [zeros(180,1), r*cosd(deg), r*sind(deg)];
tri_idx = dsearchn(mp_a,pt);
radpts = zeros(180,3);
for deg = 1:180
    radpts(deg,:) = triangle_intersect(p_a(fid.hsp.e(tri_idx(deg),:),:),pt(deg,:)) + ori;
    radpts(deg,:) = radpts(deg,:) - planenorm.*(dot(planenorm,radpts(deg,:))+d)./sum(planenorm.^2) - ori;
end
tri_idx = dsearchn(mp_a,radpts);
for deg = 1:180
    radpts(deg,:) = triangle_intersect(p_a(fid.hsp.e(tri_idx(deg),:),:),radpts(deg,:));
end
tri_idx = dsearchn(mp_a,radpts);
for deg = 1:180
    radpts(deg,:) = triangle_intersect(p_a(fid.hsp.e(tri_idx(deg),:),:),radpts(deg,:)) + ori;
end
%figure, plot3(radpts(:,1),radpts(:,2),radpts(:,3),'.'); hold on; trimesh(fid.hsp.e,fid.hsp.p(:,1),fid.hsp.p(:,2),fid.hsp.p(:,3)); hold off;

x7idx = dsearchn(radpts,x7);                % Index of X7
xzidx = dsearchn(radpts,xz);                % Index of Xz
x8idx = dsearchn(radpts,x8);                % Index of X8

lcirc = cumsum(sqrt(sum(diff(radpts(x7idx:xzidx,:)).^2,2)));
rcirc = cumsum(sqrt(sum(diff(radpts(xzidx:x8idx,:)).^2,2)));

refpoints = zeros(10,3);
% Left side
for i = 1:length(lcirc)
    fracpos = lcirc(i)/lcirc(end);
    if fracpos <= 0.1; refpoints(1,:) = radpts(x7idx+i,:);
    elseif fracpos <= 0.3; refpoints(2,:) = radpts(x7idx+i,:);
    elseif fracpos <= 0.5; refpoints(3,:) = radpts(x7idx+i,:);
    elseif fracpos <= 0.7; refpoints(4,:) = radpts(x7idx+i,:);
    elseif fracpos <= 0.9; refpoints(5,:) = radpts(x7idx+i,:);
    else; break; end;
end

% Right side
for i = 1:length(rcirc)
    fracpos = rcirc(i)/rcirc(end);
    if fracpos <= 0.1; refpoints(6,:) = radpts(xzidx+i,:);
    elseif fracpos <= 0.3; refpoints(7,:) = radpts(xzidx+i,:);
    elseif fracpos <= 0.5; refpoints(8,:) = radpts(xzidx+i,:);
    elseif fracpos <= 0.7; refpoints(9,:) = radpts(xzidx+i,:);
    elseif fracpos <= 0.9; refpoints(10,:) = radpts(xzidx+i,:);
    else; break; end;
end



%%
% Determines whether the line between a point and the origin intersects the
% specified triangle.
% @param tri_pnts the points of the triangle.
% @param pnt the point used for the line-plane intersection.
% @return isect the coordinates of the intersection if one exists.
% Otherwise an empty matrix is returned.
%%
function isect = triangle_intersect(tri_pnts, pnt)

tuv = inv([pnt' tri_pnts(2,:)'-tri_pnts(1,:)' tri_pnts(3,:)'-tri_pnts(1,:)']) * (pnt'-tri_pnts(1,:)');
%if tuv>= 0 & tuv<=1 & tuv(2)+tuv(3)<=1
    isect = pnt*(1-tuv(1));
%else
%    isect = [];
%end


%------------
function [Y,XYZ] = spm_read_vols(V,mask)
% Read in entire image volumes
% FORMAT [Y,XYZ] = spm_read_vols(V,mask)
% V    - vector of mapped image volumes to read in (from spm_vol)
% mask - implicit zero mask?
% Y    - 4D matrix of image data, fourth dimension indexes images
% XYZ  - 3xn matrix of XYZ locations returned
%_______________________________________________________________________
%
% For image data types without a representation of NaN (see spm_type),
% implicit zero masking can be used. If mask is set, then zeros are
% treated as masked, and returned as NaN.
%_______________________________________________________________________
% @(#)spm_read_vols.m	2.5 Andrew Holmes 00/03/21
% 
% Modified for NUTEEG - Daniel D.E. Wong


%-Argument checks
%-----------------------------------------------------------------------
if nargin<2, mask = 0; end
if nargin<1, error('insufficient arguments'), end

%-Image dimension, orientation and voxel size checks
%-----------------------------------------------------------------------
if length(V)>1 & any(any(diff(cat(1,V.dim),1,1),1)&[1,1,1,0])
	error('images don''t all have the same dimensions'), end
if any(any(any(diff(cat(3,V.mat),1,3),3)))
	error('images don''t all have same orientation & voxel size'), end

%-Read in image data
%-----------------------------------------------------------------------
n  = prod(size(V));			%-#images
Y = zeros([V(1).dim(1:3),n]);		%-image data matrix

for i=1:n, for p=1:V(1).dim(3)
	Y(:,:,p,i) = spm_slice_vol(V(i),spm_matrix([0 0 p]),V(i).dim(1:2),0);
end, end

%-Apply implicit zero mask for image datatypes without a NaNrep
%-----------------------------------------------------------------------
if mask
	%-Work out images without NaNrep
	im = logical(zeros(n,1));
	for i=1:n, im(i)=~spm_type(V(i).dim(4),'NaNrep'); end
	%-Mask
	Y(Y(:,:,:,im)==0)=NaN;
end

%-Return as 3D matrix if single image
%-----------------------------------------------------------------------
% if n==1; Y=Y(:,:,:,1); end

%-Compute XYZ co-ordinates (if required)
%-----------------------------------------------------------------------
if nargout>1
	[R,C,P]=ndgrid(1:V(1).dim(1),1:V(1).dim(2),1:V(1).dim(3));
    R = single(R); C = single(C); P = single(P);
	RCP = [R(:)';C(:)';P(:)'];
	RCP(4,:)=1;
	XYZ = V(1).mat*RCP;
	XYZ=XYZ(1:3,:);
end
