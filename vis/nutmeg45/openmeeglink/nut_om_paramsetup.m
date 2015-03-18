function nut_om_paramsetup(coreg,datatype,subjname,refCoord)
% This version outputs SPM-convention MRI coordinates (mm)
%%% NOTES
% FSL and BrainVisa meshes are output in mm (?)
% OpenMEEG expects meters
% additionally, MEG sensors from nutmeg need to be converted to BrainVisa
% MRI coordinates (mm), then to meters
% "refcoord" are the coordinates of the reference electrode (for EEG)

segprogram = 'brainvisa';

scalefactor = 1; % might need to be 1000 for mm <-> m conversions... but BrainVisa meshes seem to be in mm currently!

% voichunksize = 10000; % 10000 voxel chunks need ~1.4 GB/chunk
voichunksize = 20000; % 20000 voxel chunks need ~2.6 GB/chunk

global nuts st

% ctf=readCTFds('behcl_DISCR_20090427_01-epoched-std.ds');
% nuts.meg.sensorCoord = [ctf.res4.senres.pos]'
% nuts.meg.sensorOrient = [ctf.res4.senres.ori]'

sensorCoord_nuts = nuts.meg.sensorCoord;

% sensorCoord = nut_mm2voxels(nut_meg2mri(sensorCoord_nuts))/scalefactor;
sensorCoord = nut_meg2mri(sensorCoord_nuts)/scalefactor;

if(~isfield(nuts,'voxels'))
    nuts.voxels = nut_makeMNIvoi(5,'wholebrain');
end

voxels = nut_meg2mri(nuts.voxels);



%%% Fiducial file not necessary since we're computing our own transforms...
% fids = nuts.coreg.fiducials_mri_mm/1000;
% save NWfids.txt fids -ascii

% open text file for writing

switch(datatype)
    case 'MEG'
        fid = fopen([subjname '_MEGcoords.txt'],'wt');
        %% transforming orientation is ugly, but works
        sensorOrient_nuts = nuts.meg.sensorOrient;

        mat = st.vols{1}.mat;
        R=spm_imatrix(st.vols{1}.mat);
        R = spm_matrix([0 0 0 R(4:6)]);
        R = R(1:3,1:3);
        dim = st.vols{1}.dim*inv(R);
        dim_mm = abs(dim .* diag(st.vols{1}.mat(1:3,1:3)*inv(R))');
        voxsize_mm = abs(diag(st.vols{1}.mat(1:3,1:3)*inv(R))');

        % corresponding orientation transform needs to ignore translations
        sensorOrient=sensorOrient_nuts*nuts.coreg.meg2mri_tfm(1:3,1:3)';
        %%
        
        for ii=1:size(sensorCoord,1)
            fprintf(fid,'%s\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\n',[nuts.meg.rawsensor_labels{ii}],sensorCoord(ii,:),sensorOrient(ii,:));
        end
    case 'Neurospin'
        fid = fopen([subjname '_MEGcoords.txt'],'wt');
        %% transforming orientation is ugly, but works
        sensorOrient_nuts = nuts.meg.sensorOrient;

        mat = st.vols{1}.mat;
        R=spm_imatrix(st.vols{1}.mat);
        R = spm_matrix([0 0 0 R(4:6)]);
        R = R(1:3,1:3);
        dim = st.vols{1}.dim*inv(R);
        dim_mm = abs(dim .* diag(st.vols{1}.mat(1:3,1:3)*inv(R))');
        voxsize_mm = abs(diag(st.vols{1}.mat(1:3,1:3)*inv(R))');

        % corresponding orientation transform needs to ignore translations
        sensorOrient=sensorOrient_nuts*(nuts.coreg.meg2mri_tfm(1:3,1:3)';
        %%

        % hack for Neurospin data, since our test dataset doesn't have real channel labels
        for ii=1:size(sensorCoord,1)
            fprintf(fid,'%s\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\n',['nmag' num2str(ii)],sensorCoord(ii,:),sensorOrient(ii,:));
        end
    case 'EEG'
        fid = fopen([subjname '_EEGcoords.txt'],'wt');
        % add reference electrode coordinates to list
        if(exist('refCoord','var'))
            sensorCoord(end+1,:)=nut_meg2mri(refCoord)/scalefactor;
        end
        
        for ii=1:size(sensorCoord,1)
            fprintf(fid,'%s\t%.15f\t%.15f\t%.15f\n', nuts.meg.rawsensor_labels{ii}, sensorCoord(ii,:));
        end
end


fclose(fid);





% TODO: ensure that voxels are in MRI coordinates!!!
%voxels = nut_mm2voxels(nut_meg2mri(nuts.voxels)) / scalefactor;
% voxels = nut_meg2mri(nuts.voxels) / scalefactor;
numfiles = ceil(size(voxels,1)/voichunksize);

% handle dipole file
dipfile = [subjname '_voxels'];
for ii=1:numfiles
    select = ((ii-1)*voichunksize + 1) : (min((ii)*voichunksize,size(voxels,1)));
    ndip = length(select);
    % save pos with each 3D orientation
    writevoxels = [kron(voxels(select,:),ones(3,1)) , kron(ones(ndip,1),eye(3))];
    om_save_full(writevoxels,[dipfile num2str(ii) '.bin'],'binary');
end


if(0) % special voxel file construction for electrocortical stimulation project
    dipfile = 'om_ecs_pos';

    if(0) % just to remember how to load up info...
        load elecs_finalcoords
        coords = cat(1,probecoords{:});
        
        load stimdata_37864
    end

    stimselect = find(stimelecs(:,1)~=0); % remove "blank" stim trials
    stimpair_jp = stimelecs(stimselect,:);
    % coords_bv = nut_voxels2bv(nut_mm2voxels(coords));


    
    % calculate orientation of stimulating bipolar pair, and normalize
    ecs_vec = coords(stimpair_jp(:,1),:)-coords(stimpair_jp(:,2),:);
    for ii=1:size(ecs_vec,1)
        ecs_vec(ii,:) = ecs_vec(ii,:)./norm(ecs_vec(ii,:));
        
        % set dipole origin to be midpoint of stimulating bipolar pair
        ecs_pos(ii,:) = mean(coords(stimpair_jp(ii,:),:));
    end    
    
    % construct OM dipole file
    dips = [ecs_pos ecs_vec];
    om_save_full(dips,[dipfile '.bin'],'binary');
end

switch(segprogram)
    case 'brainvisa'
        nut_bvdec([subjname '_BEM_head.mesh'],'spm');
        nut_bvdec([subjname '_BEM_skull.mesh'],'spm');
        nut_bvdec([subjname '_BEM_brain.mesh'],'spm');

        nut_om_write_geom_bv(subjname);

    case 'fsl'
        if(~exist([subjname '_inskull_mesh.off'],'file'))
            disp('Executing BET...')
            if(isfield(coreg,'mriT2path'))
                system(['bet ' coreg.mripath ' ' subjname ' -A2 ' coreg.mriT2path ' -f 0.5 -g 0' -e]);
                nut_fsldec([subjname '_brain_mesh.vtk']);
            else
                system(['bet ' coreg.mripath ' ' subjname ' -A -f 0.5 -g 0']);
            end
            disp('BET done.');
        end

        nut_fsldec([subjname '_inskull_mesh.off']);
        nut_fsldec([subjname '_outskull_mesh.off']);
        nut_fsldec([subjname '_outskin_mesh.off']);
        
        nut_om_write_geom_fsl(subjname);
        warning('should write cond file here...');

    otherwise
        warning('NOT decimating meshes...');
end

system(['/usr/local/bin/om_check_geom -g ' subjname '.geom']);


function nut_om_write_geom_fsl(subjname)
fid=fopen([subjname '.geom'],'w');
fprintf(fid,'# Domain Description 1.0');
fprintf(fid,'\n');
fprintf(fid,'Interfaces 3 Mesh')
fprintf(fid,'\n');
fprintf(fid,[subjname '_inskull_mesh_small.tri']);
fprintf(fid,'\n');
fprintf(fid,[subjname '_outskull_mesh_small.tri']);
fprintf(fid,'\n');
fprintf(fid,[subjname '_outskin_mesh_small.tri']);
fprintf(fid,'\n\n');
fprintf(fid,'Domains 4')
fprintf(fid,'\n');
fprintf(fid,'Domain Scalp 1 -3');
fprintf(fid,'\n');
fprintf(fid,'Domain Brain -2');
fprintf(fid,'\n');
fprintf(fid,'Domain Air 3');
fprintf(fid,'\n');
fprintf(fid,'Domain Skull 2 -1');
fprintf(fid,'\n');
fclose(fid);


function nut_om_write_geom_bv(subjname)
fid=fopen([subjname '.geom'],'w');
fprintf(fid,'# Domain Description 1.0')
fprintf(fid,'\n');
fprintf(fid,'Interfaces 3 Mesh');
fprintf(fid,'\n');
fprintf(fid,[subjname '_BEM_brain_spm_small.mesh']);
fprintf(fid,'\n');
fprintf(fid,[subjname '_BEM_skull_spm_small.mesh']);
fprintf(fid,'\n');
fprintf(fid,[subjname '_BEM_head_spm_small.mesh']);
fprintf(fid,'\n\n');
fprintf(fid,'Domains 4')
fprintf(fid,'\n');
fprintf(fid,'Domain Scalp 1 -3');
fprintf(fid,'\n');
fprintf(fid,'Domain Brain -2');
fprintf(fid,'\n');
fprintf(fid,'Domain Air 3');
fprintf(fid,'\n');
fprintf(fid,'Domain Skull 2 -1');
fprintf(fid,'\n');
fclose(fid);


