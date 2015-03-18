function xyz_o = nut_tal2mni(xyz_i)
% loads MNI database specified in nut_defaults

global ndefaults
switch(ndefaults.mni2tal)
    case 'icbm'
        xyz_o = tal2icbm_spm(xyz_i);
    case 'brett'
        xyz_o = tal2mni(xyz_i);
    otherwise
        error('unknown MNI database specified in nut_defaults')
end
