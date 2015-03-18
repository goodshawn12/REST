function nut_bvdec(bvmeshfile,coordsys)
% decimates BrainVisa meshes, because BrainVisa's
% AimsMeshDecimation is too anal to do it properly
% input: meshfile.mesh (or .tri)
% output: meshfile_small.mesh (or .tri)
% argument: coordsys
%            'bv'  no conversion from BV convention
%            'spm' additionally transforms BV mesh to SPM-style coordinates (mm)
%            'ctf' additionally transforms BV mesh to CTF-style MEG head coords (mm)
%            'spm' and 'ctf' require Nutmeg to be open with MRI loaded/coregistered
% requires: iso2mesh toolbox (http://iso2mesh.sf.net)
%           BrainVisa matlab function (for mesh format): brainvisa/matlab/mesh/loadmesh
%           OpenMEEG matlab functions (for tri format): read_om_tri

[pathname,filename,meshfmt]=fileparts(bvmeshfile);

switch(meshfmt)
    case '.mesh'
        [bvmesh.vertices,bvmesh.faces]=loadmesh(bvmeshfile);
        bvmesh.faces = bvmesh.faces + 1;  % nerds count from 0; matlab counts from 1
    case '.tri'
        [bvmesh.vertices,bvmesh.faces]=read_om_tri(bvmeshfile);
end

bvmesh.faces = bvmesh.faces(:,[1 2 3]); % flip face order for OpenMEEG

        

% decimate to 20% of original size
[bvmeshdec.vertices,bvmeshdec.faces]=meshresample(bvmesh.vertices,bvmesh.faces,.2);

switch(coordsys)
    % WARNING: this might go funky if MRI is not isotropic...
    case 'ctf'
        bvmeshdec.vertices = nut_mri2meg(nut_voxels2mm(nut_bv2voxels(bvmeshdec.vertices)));
        filename = [filename '_ctf'];
    case 'spm'
        bvtfm = nut_read_bvtfm(bvmeshfile);
        bvmeshdec.vertices = nut_coordtfm(bvmeshdec.vertices,bvtfm);
        filename = [filename '_spm'];
end

bvmeshdec.normals = -om_normals(bvmeshdec.vertices,bvmeshdec.faces);

switch(meshfmt)
    case '.mesh'
        bvmeshdec.faces = bvmeshdec.faces - 1; %brainvisa counts from 0
        savemesh(fullfile(pathname,[filename '_small.mesh']),bvmeshdec.vertices,bvmeshdec.faces,bvmeshdec.normals);
    case '.tri'
        om_save_tri(fullfile(pathname,[filename '_small.tri']),bvmeshdec.vertices,bvmeshdec.faces);
end
