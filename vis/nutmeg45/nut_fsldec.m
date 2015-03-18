function nut_fsldec(fslmeshfile,coordsys)
% decimates FSL meshes
% input: meshfile.off
% output: meshfile_small.tri
% argument: coordsys
%            'bv'  no conversion from BV convention
%            'spm' additionally transforms BV mesh to SPM-style coordinates (mm)
%            'ctf' additionally transforms BV mesh to CTF-style MEG head coords (mm)
%            'spm' and 'ctf' require Nutmeg to be open with MRI loaded/coregistered
% requires: iso2mesh toolbox (http://iso2mesh.sf.net)
%           BrainVisa matlab function (for mesh format): brainvisa/matlab/mesh/loadmesh
%           OpenMEEG matlab functions (for tri format): read_om_tri

global st

[pathname,filename,meshfmt]=fileparts(fslmeshfile);

switch(meshfmt)
    case '.vtk'
        [fslmesh.vertices,fslmesh.faces]=read_vtk(fslmeshfile);
    case '.off'
        [fslmesh.vertices,fslmesh.faces]=readoff(fslmeshfile);
end

% fslmesh.vertices = nut_voxels2mm(fslmesh.vertices);

% FSL mesh in MRI mm, but needs to be translated to SPM origin
tfm = round(st.vols{1}.mat(1:3,1:3));
tfm(1:4,4) = st.vols{1}.mat(:,4);
fslmesh.vertices = nut_coordtfm(fslmesh.vertices,tfm);

% this will give us coord/face order (in case voxels are oriented
% differently from real mm)
tfm(:,4)=[0 0 0 1];
faceorient = nut_coordtfm([1 2 3],tfm);
if(prod(faceorient)<0)
    faceorient = abs(fliplr(faceorient));
end

faceorient = [1 3 2];

fslmesh.faces = fslmesh.faces(:,faceorient); % change face order to make OpenMEEG happy

% decimate to 20% of original size
[fslmeshdec.vertices,fslmeshdec.faces]=meshresample(fslmesh.vertices,fslmesh.faces,.7);

% fslmeshdec.normals = -om_normals(fslmeshdec.vertices,fslmeshdec.faces);

switch(meshfmt)
%     case '.mesh'
%         fslmeshdec.faces = fslmeshdec.faces - 1; %brainvisa counts from 0
%         savemesh(fullfile(pathname,[filename '_small.mesh']),fslmeshdec.vertices,fslmeshdec.faces,fslmeshdec.normals);
    case {'.off','.vtk'}
        om_save_tri(fullfile(pathname,[filename '_small.tri']),fslmeshdec.vertices,fslmeshdec.faces);
end