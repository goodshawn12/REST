

% example loading CTF data
datafile='/net/pppjz.nottingham.ac.uk/zumer/somato_demo_testing/ZumerJohanna_SEF_20031212_01.ds/ZumerJohanna_SEF_20031212_01.meg4';
sessfile='/net/pppjz.nottingham.ac.uk/zumer/somato_demo_testing/testsess.mat';

optionmode='mni';
optionmode='voxels';

switch optionmode
    case 'mni'
        mrifile='/net/cador/Anatomical_Library/JZumer/zumerN.img';
        mnifile='/net/pppjz.nottingham.ac.uk/zumer/somato_demo_testing/wzumerN.img';
        nut_create_session_file(datafile,'sessionfile',sessfile,'mrifile',mrifile,'mnifile',mnifile);
    case 'voxels'
        mrifile='/net/cador/Anatomical_Library/JZumer/zumerN.img';
        voxfile='/net/pppjz.nottingham.ac.uk/zumer/somato_demo_testing/voxels.mat';
        nut_create_session_file(datafile,'sessionfile',sessfile,'mrifile',mrifile,'voxelsfile',voxfile);
end