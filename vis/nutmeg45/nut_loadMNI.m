function MNIdm = nut_loadMNI
% loads MNI database specified in nut_defaults

global ndefaults
switch(ndefaults.mni2tal)
    case 'icbm'
        load('MNIicbm');
    case 'brett'
        load('MNIbrett');
    otherwise
        error('unknown MNI database specified in nut_defaults')
end
