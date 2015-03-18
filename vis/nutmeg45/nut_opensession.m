function nut_opensession(savedfile,partial)
% NUT_OPENSESSION
%
% The function is used to set the button and uimenu properties 
% and allows the user to work with previously stored *.mat file.
%
%

if nargin<2, partial=false; end

figh = findobj('tag','nutmegfig');
isgui = ~isempty(figh);

if ~exist('savedfile','var')
    [nutfilename, nutpathname]=uigetfile('*.mat','Open NUTMEG session...');
    if isequal(nutfilename,0)|isequal(nutpathname,0)
        return;
    end 
    savedfile=fullfile(nutpathname,nutfilename);
end

global nuts

warning('off','MATLAB:load:variableNotFound');
load(savedfile,'filename');   %assuming that 'filename' exists in workfile from older saved session (but not new NUTS way)
warning('on','MATLAB:load:variableNotFound');
if exist('filename','var')
    nuts=load(savedfile);
    clear filename
else
    clear filename
    if partial
        warning('off','MATLAB:load:variableNotFound');
        nuts=load(savedfile,'coreg','voxels','voxelsize','Lp','voxor');
        warning('on','MATLAB:load:variableNotFound');
    else
        nuts=load(savedfile);
    end
    isoldsavestyle = isfield(nuts,'nuts');
    if isoldsavestyle
        nuts=nuts.nuts;
    end
end
if isgui
    nuts.fig = figh; % update nutmeg fig handle
elseif isfield(nuts,'fig')
    nuts=rmfield(nuts,'fig'); 
end  

if(isfield(nuts,'megfile'))
    set(findobj('Tag','nut_megfile'),'String',nuts.meg.filename);
end

if(isfield(nuts,'preprocessing'))
    if(isfield(nuts.preprocessing,'beamformertype'))
        if(nuts.preprocessing.beamformertype == 1)  % legacy support... old way had numbers instead of strings
            nuts.preprocessing.beamformertype = 'Default Beamformer';
        end
    end
end

if isgui
    nut_refresh_image;  % load image into SPM
    nut_enabler;
    nut_defaults;
end
