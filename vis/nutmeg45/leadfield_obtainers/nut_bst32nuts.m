function nuts = nut_bst32nuts(lp,nuts)
% nut_bst32nuts   imports leadfield and voxel grid created with Brainstorm 3
%   nuts = nut_nut_bst32nuts(lp)
%   nuts = nut_nut_bst32nuts(lp,nuts)
%
% nuts NUTMEG session structure.
% lp   Brainstorm head model structure (obtained by exporting Head Model to Matlab
%      variable in Brainstorm 3).

nuts.voxels = lp.GridLoc' .* 1000;
nuts.voxelsize = str2num(inputdlg('Indicate approximate voxel size','Import bst head model',1,{'[5 5 5]'}));
lpbsdim = size(lp.Gain);
switch lpbsdim(2)
    case size(nuts.voxels,1)
        nuts.Lp(:,1,:) = lp.Gain;
    case 3*size(nuts.voxels,1)
        nuts.Lp(:,1,:) = lp.Gain(:,2:3:end);
        nuts.Lp(:,2,:) = lp.Gain(:,1:3:end);
        nuts.Lp(:,3,:) = lp.Gain(:,3:3:end);
    otherwise
        error('Leadfield does not match voxel grid.');
end
% nuts.voxelorient = lp.GridOrient';

nuts.meg.lsc = unique([lp.Param(:).Center]','rows') .* 1000;