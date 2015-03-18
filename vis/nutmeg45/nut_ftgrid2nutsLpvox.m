function [Lp,voxels,voxelsize] = nut_ftgrid2nutsLpvox(grid)
% [Lp,voxels,voxelsize] = nut_ftgrid2nutsLpvox(grid)
% 
% INPUT: 
% grid: output from FieldTrip 
%      grid initially can be created from ft_prepare_leadfield
%      then leadfield computed from ft_compute_leadfield and added as
%      grid.leadfield, but see Fieldtrip documentation for details
% lsc_sensor_labels:  from nuts.meg.lsc_sensor_labels, assuming you have
%      already loaded data
%
% OUTPUT: can be added directly into nuts structure, i.e.
% nuts.Lp=Lp; nuts.voxels=voxels; nuts.voxelsize=voxelsize;

if isfield(grid,'cfg') && strcmp(grid.cfg.normalize,'yes')
    error('please load a grid NOT normalized. else will interfere with inverse methods later')
end
if isfield(grid,'leadfield')
    [nk,nd]=size(grid.leadfield{grid.inside(1)});
    Lptmp=reshape([grid.leadfield{grid.inside}],nk,nd,length(grid.inside));
    % % match_str is a fieldtrip function;
    % % sensor order from ft doesn't match nm necessarily
    % % haven't decided if this necessary or not, something not working..
    % [~, chanindx] = match_str(lsc_sensor_labels, grid.cfg.channel);
    % Lp=Lptmp(chanindx,:,:);
    Lp=Lptmp;
else
    Lp=[];
end
voxels=10*grid.pos(grid.inside,:);
% ft uses cm, nm uses mm
if isfield(grid,'xgrid')
    voxelsize=[mode(10*diff(grid.xgrid)) mode(10*diff(grid.ygrid)) mode(10*diff(grid.zgrid))];
else % this is more of a guess
    ardv=abs(round(diff(voxels)));voxelsize(1:3)=mode(ardv(ardv>0));
end
        
warning('for now, please check yourself if grid.cfg.channel matches nuts.meg.sensor_labels!!');

