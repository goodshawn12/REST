function nut_iso2surf()

global st
if ~isfield(st,'vols') | isempty(st.vols{1})
    errordlg('Please load corresponding MRI first.');
    return;
end

if(~exist('vol2surf','file'))
    errordlg('The iso2mesh library is required in the MATLAB path. It can be downloaded from http://iso2mesh.sourceforge.net')
    return
end

[filename,pathname] = uigetfile({'*.img;*.nii'},'Open Segmentation Iso-Value MRI File') ;
V = spm_vol(fullfile(pathname,filename));
[Y,XYZ] = spm_read_vols(V);

lyrvals = unique(reshape(Y,[1,numel(Y)]));
lyrvals = sort(lyrvals(lyrvals>0));

prompt = 'Enter iso-values of layers to generate (outside-in):';
title = 'Layer Selection';
lines = 1;
def{1} = num2str(lyrvals);
answer = inputdlg(prompt,title,lines,def);
if isempty(answer); answer = def; end;
lyrvals = str2num(cell2mat(answer));

prompt = 'Enter maximum number of nodes (outside-in):';
title = 'Max Nodes';
lines = 1;
def{1} = num2str([1:length(lyrvals)]*1000);
answer = inputdlg(prompt,title,lines,def);
if isempty(answer); answer = def; end;
maxnode = str2num(cell2mat(answer));

% Ref - Oostendorp et al., IEEE Trans. Biomed. Eng., 2000
if length(lyrvals) == 3
    cond = [0.2025 0.0131 0.2025];
elseif length(lyrvals) == 4
    cond = [0.2025 0.0131 1.3 0.15];
else
    cond = [];
end
prompt = 'Enter tissue conductivities (outside-in):';
title = 'Tissue Conductivities';
lines = 1;
def{1} = num2str(cond);
answer = inputdlg(prompt,title,lines,def);
if isempty(answer); answer = def; end;
vol.cond = str2num(cell2mat(answer));

%floodseed = round(nut_mm2voxels(nut_mni2mri([0 0 15])));
for ii = 1:length(lyrvals)
    Ymesh = Y>=lyrvals(ii);
    %Ymesh = imfill(Ymesh,floodseed);

    %opt.keepratio=0.1; % this option is only useful when vol2mesh uses 'simplify' method
    opt.radbound=3;    % set the target surface mesh element bounding sphere be <3 pixels in radius
    opt.maxnode=maxnode(ii);
    opt.maxsurf=1;
    [node,face]=v2s(Ymesh,1,opt,'cgalsurf');
    face = face(:,1:3);
    %face = meshreorient(node,face);
    try
        node = node./repmat(st.vols{1}.private.hdr.dime.pixdim(2:4),size(node,1),1);
    catch   % The structure of st seems to be different sometimes (SPM8 vs SPM2?)
        node = node./repmat(double(st.vols{1}.private.hdr.pixdim(2:4)),size(node,1),1);
        % The double here is necessary, because pixdim is single and will
        % convert cs.vertices to single as well. This in turn will kill matlab in
        % reducepatch!
    end
    node = nut_voxels2mm(node);
    
    vol.bnd(ii).vertices = node;
    vol.bnd(ii).faces = face;
end

vol = decouplesurf(vol);

% for ii = 1:length(lyrvals)
%     node = vol.bnd(ii).vertices;
%     face = vol.bnd(ii).faces;
%     
%     face = meshreorient(node,face);
%     
%     % Find a face with highest node z-value - we'll use this to figure out
%     % the orientations
%     [~,I] = max(node(:,3));   % Get node on top of head
%     topface = face(find(face(:,1)==I | face(:,2)==I | face(:,3)==I),:); topface = topface(1,:);
%     
%     nvec = cross(node(topface(2),:)-node(topface(1),:),node(topface(3),:)-node(topface(1),:));
%     tricenter = (node(topface(1),:)+node(topface(2),:)+node(topface(3),:))/3;
%     th = dot(nvec,tricenter);
%     if th > 0
%         face = [face(:,3) face(:,2) face(:,1)];
%     end
%     vol.bnd(ii).vertices = node;
%     vol.bnd(ii).faces = face;
% end
for ii = 1:length(lyrvals)
    [vol.bnd(ii).vertices, vol.bnd(ii).faces] = surfreorient(vol.bnd(ii).vertices, vol.bnd(ii).faces);
    vol.bnd(ii).faces = vol.bnd(ii).faces(:,[3 2 1]);
end

nut_trimesh(vol); pause(1)

[file,path] = uiputfile('*_vol.mat','Save Surface Mesh');
if ~isequal(file,0) & ~isequal(path,0)
    save(strcat(path,file),'vol');
end



%%
% Take care of intersecting surfaces.
%%
function vol = decouplesurf(vol)

for ii = 1:length(vol.bnd)-1
    % Despite what the instructions for surfboolean says, surfaces should
    % be ordered from inside-out!!
    [newnode, newelem] = surfboolean(vol.bnd(ii+1).vertices,vol.bnd(ii+1).faces,'decouple',vol.bnd(ii).vertices,vol.bnd(ii).faces);
    vol.bnd(ii+1).faces = newelem(newelem(:,4)==2,1:3) - size(vol.bnd(ii+1).vertices,1);
    vol.bnd(ii+1).vertices = newnode(newnode(:,4)==2,1:3);
end