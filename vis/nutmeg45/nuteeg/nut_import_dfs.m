%%
% Imports head layers from BrainSuite for use with EEG head model tools.
% @author Daniel D.E. Wong
%%
function nut_import_dfs

% Make sure MRI is loaded so we can convert units to mm
global st
if ~isfield(st,'vols') | isempty(st.vols{1})
    error('Please load corresponding MRI first.');
    return;
end

ns = 0; % Number of surfaces loaded
path = '';
while 1
    [file,path] = uigetfile('*.dfs','Open DFS (Outside-In)',path);
    if ~isequal(file,0) & ~isequal(path,0)
        ns = ns + 1;
        dfsfiles{ns} = strcat(path,file);
        fprintf('[%d] %s\n',ns,file);
        vol.filenames{ns} = file;
        
        prompt = ['Decimate faces down to:'];
        title = 'Import DFS';
        lines = 1;
        def{1} = num2str(1000*ns);
        answer = inputdlg(prompt,title,lines,def);
        if ~isempty(answer); d(ns) = str2num(cell2mat(answer)); else; d(ns) = 1; end;
        if d(ns) > 5120; disp('Maximum faces allowed for decimation is 5120.'); end;
    else
        break;
    end
end
if ~ns; return; end;

prompt = 'Enter layer conductivities:';
title = 'Import DFS';
lines = 1;
if ns >= 4; def{1} = '[0.33 0.0041 1.79 0.33]'; else; def{1} = '[0.33 0.0041 0.33]'; end;
answer = inputdlg(prompt,title,lines,def);
if ~isempty(answer); vol.cond = str2num(cell2mat(answer)); else, return, end;

for i = 1:ns
    fprintf('Reading %s...\n',vol.filenames{i});
    dfs = readdfs(dfsfiles{i});
    
    % Convert to MRI [mm] dimensions
    try
        dfs.vertices = dfs.vertices./repmat(st.vols{1}.private.hdr.dime.pixdim(2:4),size(dfs.vertices,1),1);
    catch   % The structure of st seems to be different sometimes (SPM8 vs SPM2?)
        dfs.vertices = dfs.vertices./repmat(double(st.vols{1}.private.hdr.pixdim(2:4)),size(dfs.vertices,1),1);
        % The double here is necessary, because pixdim is single and will
        % convert cs.vertices to single as well. This in turn will kill matlab in
        % reducepatch!
    end
    dfs.vertices = nut_voxels2mm(dfs.vertices);
    
    if(d(i) ~= 1)
        dfs = nut_shrinkfitsphere(reducepatch(dfs,d(i)),300,d(i));
    end
    vol.bnd(i).vertices = dfs.vertices;
    vol.bnd(i).faces = dfs.faces;
end

[file,path] = uiputfile('_vol.mat','Save to MAT File');
if ~isequal(file,0) & ~isequal(path,0)
    save(strcat(path,file),'vol');
end

nuts.coreg.volfile = fullfile(path,file);
