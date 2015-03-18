%%
% Displays triangular mesh, with option of showing normal vectors.
% @author Daniel D.E. Wong
%%
function nut_trimesh(vol)
global nuts
if nargin<1
    def_path = [];
    if isfield(nuts.coreg,'volfile')
        [def_path,~,~] = fileparts(nuts.coreg.volfile);
    end
    [file,path] = uigetfile('*_vol.mat','Open Mesh File',def_path);
    if ~isequal(file,0)
        load(strcat(path,file));
    else
        return;
    end
end

answer_n = questdlg('Show normals (this will take time)?','Display Mesh','Yes','No','No');
if isfield(nuts,'meg') && isfield(nuts.meg,'sensorCoord')
    answer_sl = questdlg('Show sensor locations?','Display Mesh','Yes','No','No');
else
    answer_sl = 'No';
end

nMeshes = length(vol.bnd);
figure('units','normalized','position',[0.01 0.5 nMeshes*0.24 0.35]);

for m=1:nMeshes
% if nMeshes>1
%     prompt = ['Which mesh to display? (1-' num2str(nMeshes) ')'];
%     title = 'Display Mesh';
%     lines = 1;
%     def{1} = '1';
%     answer = inputdlg(prompt,title,lines,def);
%     if ~isempty(answer); m = str2num(cell2mat(answer)); else; m = 1; end;
%     if m > nMeshes | m < 1; m = 1; end;
% else
%     m=1;
% end

% while 1
%     answer = questdlg(['Increase density? (' num2str(size(vol.bnd(m).faces,1),1) 'faces)'],'Display Mesh','Yes','No','No');
%     if strcmp(answer,'Yes')
%         vol.bnd(m) = mesh_refine_tri4(vol.bnd(m));
%     else
%         break;
%     end
% end

    %mesh = PrepareTriangleMesh(nut_mri2meg(vol.bnd(m).vertices),vol.bnd(m).faces(:,[1 3 2]));
    
    % Handle FieldTrip style vol
    if isfield(vol.bnd(m),'tri'); vol.bnd(m).faces=vol.bnd(m).tri; end;
    if isfield(vol.bnd(m),'pnt'); vol.bnd(m).vertices=vol.bnd(m).pnt; end;
    
    nodes = nut_mri2meg(vol.bnd(m).vertices);
    faces = vol.bnd(m).faces;

    subplot(1,nMeshes,m);
    trimesh(faces,nodes(:,1),nodes(:,2),nodes(:,3),'EdgeColor','k','FaceColor','k');
    shading faceted
    view([-209 24]); 
    axis equal
    axis off
    if strcmp(answer_sl,'Yes')
        hold on;
        plot3(nuts.meg.sensorCoord(:,1),nuts.meg.sensorCoord(:,2),nuts.meg.sensorCoord(:,3),'bo','MarkerFaceColor','b');
    end
    if strcmp(answer_n,'Yes')
        %nn = nodeNormals(nodes,faces);
        [fc,fn] = faceNormals(nodes,faces);
        for i = 1:size(fc,1)
            hold on;
            vectarrow(fc(i,:),fc(i,:)+10*fn(i,:));
        end
    end
    hold off;
    rotate3d on
end
colormap(pink);

%--------------------
function [fc,fn] = faceNormals(nodes,faces)

% Find face normals
p1 = nodes(faces(:,1),:);
p2 = nodes(faces(:,2),:);
p3 = nodes(faces(:,3),:);
vec12 = p2-p1;
vec13 = p3-p1;
fn = cross(vec13,vec12,2);
fn = fn./repmat(nut_rownorm(fn),1,3);
fc = (p1+p2+p3)/3;

% Find node normals
% nn = zeros(size(nodes,1),3);
% for ii = 1:size(faces,1)
%     nn(faces(ii,:),:) = nn(faces(ii,:),:) + fn(faces(ii,:),:);
% end
% nn = nn./repmat(nut_rownorm(nn),1,3);

