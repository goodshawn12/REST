%%
% Imports a cortical surface from BrainSuite for use with cortical surface
% constraint/projection algorithms.
% @author Daniel D.E. Wong
%%
function nut_import_corticalsurface_dfs

% Make sure MRI is loaded so we can convert units to mm
global st
if ~isfield(st,'vols') | isempty(st.vols{1})
    error('Please load corresponding MRI first.');
    return;
end

[file,path] = uigetfile('*.dfs','Open DFS (Outside-In)');
if isequal(file,0) | isequal(path,0); return; end;

cs = readdfs(strcat(path,file));
% Convert to MRI [mm] dimensions
try
    cs.vertices = cs.vertices./repmat(st.vols{1}.private.hdr.dime.pixdim(2:4),size(cs.vertices,1),1);
catch   % The structure of st seems to be different sometimes (SPM8 vs SPM2?)
    cs.vertices = cs.vertices./repmat(double(st.vols{1}.private.hdr.pixdim(2:4)),size(cs.vertices,1),1);
    % The double here is necessary, because pixdim is single and will
    % convert cs.vertices to single as well. This in turn will kill matlab in
    % reducepatch!
end
cs.vertices = nut_voxels2mm(cs.vertices);

prompt = ['Decimate faces down to:'];
title = 'Import DFS';
lines = 1;
def{1} = num2str(size(cs.faces,1));
answer = inputdlg(prompt,title,lines,def);
if ~isempty(answer); d = str2num(cell2mat(answer)); else; d = 1; end;

if(d ~= 1)
    cs = reducepatch(cs,d);
end

[file,path] = uiputfile('_cs.mat','Save to MAT File');
if ~isequal(file,0) & ~isequal(path,0)
    save(strcat(path,file),'cs');
end