% Computes lead field/potential using OpenMEEG.
% To recalculate with new sensor positions, delete:
%  *._sensorcoords.txt
%  *._h2mm.bin
%  *._ds2mm.bin
%  *._meeggain*.bin
% If files were moved to another directory, delete *.geom.
function [Lfp,voxels]=nut_om_compute_lfp(voxelsize)
global nuts ndefaults

CPU_LIM = feature('numCores');
OPENMEEG_PATH = '/usr/local/bin/';
VOXCHUNKSIZE = 20000;

LD_LIBRARY_PATH='';
if ispc
    warning('Sorry, Windows is not yet supported');
elseif isunix
    LD_LIBRARY_PATH = 'export LD_LIBRARY_PATH=/usr/lib:/usr/local/lib; ';   % Matlab screws with the LD_LIBRARY_PATH variable
end

CPU_LIM = ['export OMP_NUM_THREADS=', num2str(CPU_LIM), '; '];


if isfield(nuts.coreg,'volfile') && exist(nuts.coreg.volfile,'file')
    load(nuts.coreg.volfile);
    [path,file,ext] = fileparts(nuts.coreg.volfile); %path(end+1) = '/';
%     [path,file,ext] = fileparts(nuts.coreg.volfile);
%     [file,path] = uigetfile('*_vol.mat','Open Head Surface Mesh File',path);
%     load(strcat(path,file));
%     [path,file,ext] = fileparts([path file]); path(end+1) = '/';
else
    [file,path] = uigetfile('*_vol.mat','Open Head Surface Mesh File');
    load(fullfile(path,file));
    [path,file,ext] = fileparts([path file]); %path(end+1) = '/';
end

%file=regexpi(nuts.coreg.volfile,'.mat','split'); file=file{1}
%[path,file,ext] = fileparts(nuts.coreg.volfile); %path(end+1) = '/';
basefile=regexpi(file,'_vol','split'); basefile=basefile{1};

sensorFile = fullfile(path, [basefile '_sensorcoords.txt']);
if exist(sensorFile,'file')
    disp('Sensor coordinate file already exists. Skipping...')
else
    disp('Writing sensor coordinates...')
    fid = fopen(sensorFile,'w');
    sensorCoord = nut_meg2mri(nuts.meg.sensorCoord);
    if isfield(nuts.meg,'eegflag') & nuts.meg.eegflag
        % add reference electrode coordinates to list
        if isfield(nuts.meg,'refCoord')
            sensorCoord(end+1,:)=nut_meg2mri(refCoord);
        end

        for ii=1:size(nuts.meg.sensorCoord,1)
            fprintf(fid,'%s\t%.15f\t%.15f\t%.15f\n', nuts.meg.sensor_labels{ii}, sensorCoord(ii,:)/1000);
        end
    else    % MEG
        % Convert to MRI coords
        A = nuts.coreg.meg2mri_tfm;
        A(1:3,4)=0;                     % Get rid of translations
        sensorOrient = nut_coordtfm(nuts.meg.sensorOrient,A);
        sensorOrient = sensorOrient./repmat(nut_rownorm(sensorOrient),1,3); % Just in case

        for ii=1:size(sensorCoord,1)
            fprintf(fid,'%s\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\n',nuts.meg.rawsensor_labels{ii},sensorCoord(ii,:)/1000,sensorOrient(ii,:));
        end
    end
    fclose(fid);
end

condFile = fullfile(path, [basefile '.cond']);
if exist(condFile,'file')
    disp('Conductivity file already exists. Skipping...')
else
    disp('Writing conductivity file...')
    write_cond(vol,condFile);
end

geomFile = fullfile(path, [basefile '.geom']);
if exist(geomFile,'file')
    disp('Geometry descriptor file already exists. Skipping...')
else
    disp('Writing geometry descriptor file...')
    write_geom(vol,geomFile,fullfile(path, basefile));
end

disp('Writing OpenMEEG mesh files...')
% Only need to compute normals for .mesh files, but we're using .tri files
% for ii = 1:length(vol.bnd)  % Mesh is already in MRI coords
%     % Compute triangle normals
%     vol.bnd(ii).normals = zeros(size(vol.bnd(ii).faces));
%     v12 = vol.bnd(ii).vertices(vol.bnd(ii).faces(:,2),:) - vol.bnd(ii).vertices(vol.bnd(ii).faces(:,1),:);
%     v13 = vol.bnd(ii).vertices(vol.bnd(ii).faces(:,3),:) - vol.bnd(ii).vertices(vol.bnd(ii).faces(:,1),:);
%     facenormals = cross(v12,v13);
%     facenormals = facenormals./repmat(sqrt(sum(facenormals.^2,2)),1,3);    % Compute unit normal
%     
%     % Compute node normals
%     nodenormals = zeros(size(vol.bnd(ii).vertices,1),3);
%     for jj = 1:size(vol.bnd(ii).faces,1)
%         nodenormals(vol.bnd(ii).faces(jj,1),:) = nodenormals(vol.bnd(ii).faces(jj,1),:) + facenormals(jj,:);
%         nodenormals(vol.bnd(ii).faces(jj,2),:) = nodenormals(vol.bnd(ii).faces(jj,2),:) + facenormals(jj,:);
%         nodenormals(vol.bnd(ii).faces(jj,3),:) = nodenormals(vol.bnd(ii).faces(jj,3),:) + facenormals(jj,:);
%     end
%     vol.bnd(ii).normals = nodenormals./repmat(sqrt(sum(nodenormals.^2,2)),1,3);    % Compute unit normal
% end
write_mesh(vol,path,basefile);

disp('Validating mesh...')
[stat res] = system([LD_LIBRARY_PATH, CPU_LIM, OPENMEEG_PATH, 'om_check_geom -g ', geomFile]);
if stat > 0
    warning([res, '...quitting']);
    return
end

disp('Writing dipole file...')
voxels=nut_meg2mri(nut_make_voxels(voxelsize));
chunks = ceil(size(voxels,1)/VOXCHUNKSIZE);
dipFile = cell(chunks,1);
for ii = 1:chunks
    dipFile{ii} = fullfile(path, [basefile '_voxels' num2str(ii) '.bin']);
    if exist(dipFile{ii},'file')
        fprintf('\t%s already exists. Skipping...\n', dipFile{ii});
        break;
    else
        voxidx = ((ii-1)*VOXCHUNKSIZE + 1) : (min((ii)*VOXCHUNKSIZE,size(voxels,1)));
        writevoxels = [kron(voxels(voxidx,:),ones(3,1)) , kron(ones(length(voxidx),1),eye(3))];
        om_save_full(writevoxels/1000,dipFile{ii},'binary');
    end
end


disp('--------------------------------------')

hmFile = fullfile(path, [basefile '_hm.bin']);
if exist(hmFile,'file')
    disp('Head matrix already exists. Skipping...')
else
    disp('Building head matrix')
    system([LD_LIBRARY_PATH, CPU_LIM , fullfile(OPENMEEG_PATH, 'om_assemble'), ' -hm ', geomFile, ' ', condFile, ' ', hmFile]);
end

if isfield(nuts.meg,'eegflag') & nuts.meg.eegflag
    ohmicFile = fullfile(path, [basefile '_h2em.bin']);
    cmd = '-h2em';
else
    ohmicFile = fullfile(path, [basefile '_h2mm.bin']);
    cmd = '-h2mm';
end
if exist(ohmicFile,'file')
    disp('Ohmic current file already exists. Skipping...')
else
    disp('Calculating Contribution of Ohmic Currents')
    system([LD_LIBRARY_PATH, CPU_LIM , fullfile(OPENMEEG_PATH, 'om_assemble'), ' ', cmd, ' ', geomFile, ' ', condFile, ' ' , sensorFile, ' ' , ohmicFile]);
end

if ~isfield(nuts.meg,'eegflag') || ~nuts.meg.eegflag
    disp('Contribution of all sources to the MEG sensors')
    scFile = cell(chunks,1);
    for ii = 1:chunks
        scFile{ii} = fullfile(path, [basefile '_ds2mm' num2str(ii) '.bin']);
        if exist(scFile{ii},'file')
            fprintf('\t%s already exists. Skipping...\n',scFile{ii})
            break;
        else
            system([LD_LIBRARY_PATH, CPU_LIM , fullfile(OPENMEEG_PATH, 'om_assemble'), ' -ds2mm ', dipFile{ii} ,' ', sensorFile, ' ' , scFile{ii}]);
        end
    end
end

disp('Putting it all together. Creating the BEM might take a while')
bemFile = cell(chunks,1);
for ii = 1:chunks
    if isfield(nuts.meg,'eegflag') & nuts.meg.eegflag
        bemFile{ii} = fullfile(path, [basefile '_eeggain' num2str(ii) '.bin']);
    else
        bemFile{ii} = fullfile(path, [basefile '_meggain' num2str(ii) '.bin']);
    end
    
    if exist(bemFile{ii},'file')
        fprintf('/t%s already exists. Skipping...\n', bemFile{ii});
        break;
    elseif isfield(nuts.meg,'eegflag') & nuts.meg.eegflag
        system([LD_LIBRARY_PATH, CPU_LIM , fullfile(OPENMEEG_PATH, 'om_gain'), ' -EEGadjoint ', geomFile, ' ', condFile, ' ', dipFile{ii},' ', hmFile, ' ', ohmicFile, ' ', bemFile{ii}]);
    else
        system([LD_LIBRARY_PATH, CPU_LIM , fullfile(OPENMEEG_PATH, 'om_gain'), ' -MEGadjoint ', geomFile, ' ', condFile, ' ', dipFile{ii},' ', hmFile, ' ', ohmicFile, ' ', scFile{ii},' ', bemFile{ii}]);
    end
end

% Import lead field/potential
if isfield(nuts.meg,'eegflag') & nuts.meg.eegflag; eegflag = 1; else eegflag = 0; end;
[g, voxels] = import_gain(path, basefile, eegflag);
voxels = nut_mri2meg(voxels*1000);
if isfield(nuts.meg,'chanmixMtx'); chanmixMtx = nuts.meg.chanmixMtx{1}; else; chanmixMtx=[]; end;
Lfp = computeLfp(g, size(nuts.meg.sensorCoord,1),chanmixMtx);


function write_cond(vol,filename)
fid=fopen(filename,'w');
fprintf(fid,'# Properties Description 1.0 (Conductivities)\n');
tissues = {'Scalp\t%f\n'; 'Skull\t%f\n'; 'CSF\t%f\n'; 'Brain\t%f\n'};
if length(vol.cond)==3; tissues = tissues([1 2 4]); end;
fprintf(fid,'Air\t0\n');
for ii=1:length(vol.cond)
    fprintf(fid,tissues{ii},vol.cond(ii));
end
fclose(fid);


function write_geom(vol,filename,basepathfile)
fid=fopen(filename,'w');
fprintf(fid,'# Domain Description 1.0\n');
fprintf(fid,'Interfaces %i Mesh\n',length(vol.cond));
tissues={'_scalp.tri\n'; '_skull.tri\n'; '_csf.tri\n'; '_brain.tri\n'};
if length(vol.cond)==3; tissues = tissues([1 2 4]); end;
for ii = 1:length(vol.cond)
    fprintf(fid,[basepathfile tissues{ii}]);
end
fprintf('\n');
fprintf(fid,'Domains %i\n',length(vol.cond)+1);
domains={'Scalp'; 'Skull'; 'CSF'; 'Brain'};
if length(vol.cond)==3; domains = domains([1 2 4]); end;
fprintf(fid,'Domain Air %i\n',1);
for ii = 1:length(vol.cond)
    if ii < length(vol.cond)
        fprintf(fid,['Domain ' domains{ii} ' %i -%i\n'],ii+1,ii);
    else
        fprintf(fid,['Domain ' domains{ii} ' -%i\n'],ii);
    end
end
fclose(fid);


function write_mesh(vol,path,basefile)
tissues={'_scalp'; '_skull'; '_csf'; '_brain'};
if length(vol.cond)==3; tissues = tissues([1 2 4]); end;
for ii = 1:length(vol.cond)
    meshFile = fullfile(path, [basefile tissues{ii} '.tri']);
    if exist(meshFile,'file')
        fprintf('\t%s already exists. Skipping...\n', meshFile);
        break;
    else
        om_save_tri(fullfile(path, [basefile tissues{ii} '.tri']),vol.bnd(ii).vertices/1000,vol.bnd(ii).faces);
        %savemesh([path basefile tissues{ii} '.mesh'],vol.bnd(ii).vertices/1000,vol.bnd(ii).faces-1,-vol.bnd(ii).normals);
    end
end


function Lfp = computeLfp(g, nsensors, chanmixMtx)
if isempty(chanmixMtx)
    chanmixMtx=eye(nsensors);
    if nsensors < size(g,1)
        % A reference sensor is included in g - rereference
        chanmixMtx(:,end+1) = -1;
    end
end
Lfp = chanmixMtx*g;
Lfp = reshape(Lfp,size(Lfp,1),3,size(Lfp,2)/3); % Mchannels x 3 orientations x Nvoxels


function [g, voxels] = import_gain(path, basefile, eegflag)
if eegflag
    omgainfiles = dir(fullfile(path, [basefile '_eeggain*.bin']));
else
    omgainfiles = dir(fullfile(path, [basefile '_meggain*.bin']));
end
omvoxfiles = dir(fullfile(path, [basefile '_voxels*.bin']));

g=[];
voxels=[];

% join gain/voxel files
% [openmeeg calculation may have been split for memory reasons]
for ii=1:length(omgainfiles)
    g = [g om_load_full(fullfile(path, omgainfiles(ii).name),'binary')];
    voxels = [voxels;om_load_full(fullfile(path, omvoxfiles(ii).name),'binary')];
end
voxels = voxels(1:3:end,1:3);