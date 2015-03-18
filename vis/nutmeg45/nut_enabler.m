function nut_enabler(handles)
% figures out what information is loaded into NUTMEG and
% enables/disables GUI buttons/menus as necessary

global nuts
if(~exist('handles','var'))
    if ~isempty(gcbf)  % will be empty when called from command line
        handles = guihandles(findobj('tag','nutmegfig'));
    else
        handles = guihandles(nuts.fig);
    end
end

% get handlenames of all buttons and menus, place in enableditems
% i.e., begin with everything enabled
handlenames = fieldnames(handles);
enableditems = {};
for i = 1:length(handlenames)
    % overly complicated way of testing if the handlename contains "button" or "menu"
    if(~isempty(regexpi(handlenames{i},'button')) || ~isempty(regexpi(handlenames{i},'menu')))
        enableditems{end+1} = handlenames{i};
    end
end

% now selectively disable:
disableditems = {};
if ~isfield(nuts,'meg')|~isfield(nuts.meg,'sensorCoord')  % NUTEEG mod
    disableditems = union(disableditems,{  'nut_leadfield_button',...
                                           'nut_displayMEG_button',...
                                           'nut_beamforming_button',...
                                           'nut_plotsensors_menu',...
                                           'nut_eegmralign_menu',...
                                           'nut_eeg_multisphere_menu',...
                                           'nut_export_eeg_coords_menu'});     %'nut_tfbf_button',...
else
    if ~isfield(nuts.meg,'system')
        nuts.meg.system='unknown'; %legacy support
    end
%     if ~isfield(nuts.meg,'lsc')
%         disableditems = union(disableditems,{  'nut_leadfield_button',...
%             });
%     end
    if(~isfield(nuts,'VOIvoxels') & ~isfield(nuts.coreg,'norm_mripath') & ~strcmp(nuts.meg.system,'Neuromag'))
        disableditems = union(disableditems,{  'nut_leadfield_button',...
            });
    end
end

% NUTEEG mod
if (~isfield(nuts,'meg') | ~isfield(nuts.meg,'eegflag'))
    disableditems = union(disableditems,{   'nut_import_eeg_coords_menu'  });
end

if(~isfield(nuts,'Lp') )
    disableditems = union(disableditems,{  'nut_beamforming_button',...
                                           'nut_field_viewer_menu',...
                                           'nut_tfbf_button'                }); %'nut_tfbf_button',...
end

% if(~isfield(nuts,'VOIvoxels'))
%     disableditems = union(disableditems,{  'nut_leadfield_button',...
%                                            'nut_beamforming_button'        });
% end
if(~isfield(nuts,'voxels'))
    disableditems = union(disableditems,{  'nut_beamforming_button',...
                                           'nut_tfbf_button'                });   %'nut_tfbf_button'
end


% remove explicitly disabled items from enableditems
enableditems = setdiff(enableditems,disableditems);

% enable enableditems
for i=1:length(enableditems)
    set(handles.(enableditems{i}),'Enable','On')
end

% disable disableditems
for i=1:length(disableditems)
    set(handles.(disableditems{i}),'Enable','Off')
end

% set MEG dataset name
if(isfield(nuts,'meg'))
    if(isfield(nuts.meg,'filename'))
        set(handles.nut_megfile,'String',nuts.meg.filename);
    end
end