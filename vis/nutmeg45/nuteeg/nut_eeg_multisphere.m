%%
% Finds localsphere centers for each EEG sensor using the sensor
% coordinates.
% @param sensorCoords a matrix of sensor coordinates.
% @param XYZ head surface points.
% @param Nxyz head surface point normals (pointing inward).
% @author Daniel D.E. Wong
%%
function nut_eeg_multisphere

global nuts

% Parameters
params.rinit = [30 90]; % Simplex initiallization values
params.rdefault = 100;  % Default radius if simplex search is outside rmin-rmax range
params.rmax = 150;
params.rmin = 80;

% Make sure a reference Lp solution is already loaded into nuts struct
if ~isfield(nuts,'Lp'); error('No reference lead potential solution loaded'); return; end;

% Load mesh
[file,path] = uigetfile('*_vol.mat','Open Mesh File');
if isequal(file,0); error('No Mesh file selected'); return; end;
load(strcat(path,file));
% Convert to MEG coordinates
for i = 1:length(vol.bnd)
    vol.bnd(i).vertices = nut_mri2meg(vol.bnd(i).vertices);
end
scalpmesh = PrepareTriangleMesh(vol.bnd(1).vertices,vol.bnd(1).faces(:,[1 3 2]));

h = waitbar(0,'Fitting Localspheres: 0%');
tic
for sens = 1:size(nuts.meg.sensorCoord,1)
    fprintf('Fitting local sphere for sensor: %d\n',sens);
    
    % Initialize simplex
    sensorCoord = nuts.meg.sensorCoord(sens,:);
    normvec = scalpmesh.un(dsearchn(scalpmesh.mp,sensorCoord),:);
    r0 = zeros([1,1,2]);
    r0(:,:,1) = params.rinit(1); r0(:,:,2) = params.rinit(2);
    r = nut_simplexSearch(@sphereEval,r0,[],nuts.voxels,sensorCoord,normvec,vol,squeeze(nuts.Lp(sens,:,:)));
    if r > params.rmax | r < params.rmin; r = params.rdefault; end;
    
    lsc(sens,:) = sensorCoord-normvec*r;
    fprintf('[%i %i %i] radius-to-sensor:%i\n\n',lsc(sens,1),lsc(sens,2),lsc(sens,3),r);
    
    waitbar(sens/size(nuts.meg.sensorCoord,1),h,['Fitting Localspheres: ' num2str(sens/size(nuts.meg.sensorCoord,1)*100,3) '%']);
end
toc
delete(h)
nuts.meg.lsc = lsc;
nuts.meg.lsc_sensor_labels = nuts.meg.sensor_labels;


function cost = sphereEval(r, voxels, sensorCoord, normvec, vol, LpBEM)

LpBEM = reshape(LpBEM',prod(size(LpBEM)),1);
lsc = sensorCoord-normvec*r;
[lx1,ly1,lz1] = nut_bergsphleadp(voxels,sensorCoord,lsc,vol);
LpSph = [lx1; ly1; lz1];

cost = norm(abs(LpSph)-abs(LpBEM));