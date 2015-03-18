function xyz_o = nut_mni2tal(xyz_i)
% loads MNI database specified in nut_defaults

global ndefaults
switch(ndefaults.mni2tal)
    case 'icbm'
        xyz_o = icbm_spm2tal(xyz_i);
    case 'brett'
        xyz_o = mni2tal(xyz_i);
    otherwise
        error('unknown MNI database specified in nut_defaults')
end
