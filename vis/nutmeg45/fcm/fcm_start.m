function fcm_start(varargin)
% FCM_START  sets configuration for FCM.
%
% Usages:
%   FCM_START GUI
%       opens a GUI for selection of configuration parameters.
%   FCM_START FILENAME 
%       sets the configuration saved in FILENAME.
%   FCM_START -MAKEDEFAULT FILENAME
%       sets the configuration saved in FILENAME and makes it the default.
%   FCM_START 
%       sets the default configation.
%
% See also  FCM_CONFIG_GUI

global fuse

defaultfile=[fileparts(which('fcm_start.m')) filesep 'params' filesep 'fcmConfig.ini'];
switch nargin
    case 0
        fn=textread(defaultfile,'%c','headerlines',1)';
    case 1
        if strcmpi(varargin{1},'gui')
            fcm_config_gui
            return
        else
            fn= varargin{1};
        end
    case 2
        flagidx = strmatch('-makedefault',lower(varargin));
        if isempty(flagidx), error('Invalid syntax.'), end
        fn=varargin{setdiff([1 2],flagidx)};
    otherwise
        error('Too many input arguments.')
end

load([fileparts(which('fcm_start.m')) filesep 'params' filesep fn]);

if nargin==2
    delete(defaultfile);
    fid = fopen(defaultfile,'wt');
    fprintf(fid,'[DefaultConfigFilename]\n');
    fprintf(fid,'%s\n',varargin{setdiff([1 2],flagidx)});
    fclose(fid);
    fprintf('New default file is %s.\n',varargin{setdiff([1 2],flagidx)})
end

fuse
   
    
    
