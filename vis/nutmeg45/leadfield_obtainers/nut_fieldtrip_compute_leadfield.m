function [nuts,lf_ftu]=nut_fieldtrip_compute_leadfield(nuts,sens);
% calls Fieldtrip or MJB EEG lead field code
% Fieldtrip want coords in cm, je pense.
% MJB wants coords in m, alors!
%
% does FT want 'meg' coords or 'mri'? going with 'meg' for now, as long as consistent
% likewise for mjb code
%
% need to compute nuts.meg.sphererad (in mm) beforehand

% global nuts

vol.o=nuts.meg.lsc/10; 
vol.c=[1 0.0125 1];  % conductivity ratios, mjb and fieldtrip both use these values

%get sphere radius from *.hdm or other
% nuts.meg.sphererad=92.1; %dan
% nuts.meg.sphererad=88.3; %claire
% nuts.meg.sphererad=87.0; % johanna
% nuts.meg.sphererad=87.2; % anniek
% bound=[0.87 0.92 1.00]; % radii of spheres, used by mjb
bound=[0.88 0.92 1.00]; % radii of spheres, used by FieldTrip
vol.r=nuts.meg.sphererad*bound/10; 


sens.pnt=nuts.meg.sensorCoord/10;
% sens.type='electrodes'; % for EEG
% sens.type='gradiometers'; % for MEG

if ~isfield(nuts,'voxels')
    voxels=double(nut_make_voxels(nuts.voxelsize(1)));
    nuts.voxels=voxels;
    clear voxels
end

% fieldtrip's code: have it in your path, download from fieldtrip website.
ftcode=which('ft_compute_leadfield');
if isempty(ftcode)
%     error('Please place Fieldtrip code in your path (including symbolic links to private directories)');
    error('Please place Fieldtrip main directory in your path, and call ft_defaults to set other sub-paths');
end

ft_defaults;
lf = ft_compute_leadfield(double(nuts.voxels)/10,sens,vol,'reducerank','no','normalize','no');
lf_ftu=reshape(lf,size(lf,1),3,size(nuts.voxels,1));
clear lf
% this below does same as setting cn=1 in nutmeg.
% lf = compute_leadfield(double(nuts.voxels)/10,sens,vol,'reducerank','no','normalize','column');
% lf_ftc=reshape(lf,size(nuts.meg.data,2),3,size(nuts.voxels,1));
% clear lf


%% matt's code
if(0)
    Dall = [1 0 0;0 1 0;0 0 1].*1e-9;
    for ii=1:size(nuts.voxels,1)
        for or=1:3
            D = Dall(or,:);
            lf(or,:,ii) = EEG_ScalarForward(nuts.voxels(ii,:)/1000,D,nuts.meg.sensorCoord/1000,nuts.meg.lsc/1000,nuts.meg.sphererad/1000,3,bound,vol.c);
        end
    end
    lf_mb=permute(lf,[2 1 3]);
    clear lf
end

