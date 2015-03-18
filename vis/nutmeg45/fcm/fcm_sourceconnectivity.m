function fcm_sourceconnectivity(confilename,varargin)
% FCM_SOURCECONNECTIVITY  wrapper function for calculation of source functional connectivity
%   on a single local computer.
%
% The usage depends on your configuration (which you set with fcm_start or
% fcm_config_gui).
% Call fcm_sourceconnectivity without input to see the usage for your current setting. 

global fuse

if isempty(fuse), error('You must set the configuration first (with fcm_start or fcm_config_gui)'), end

fcn = ['s' fuse.funconn fuse.datatype];
if nargin<3, 
    fprintf(helptxt(fcn))
    return 
end

feval(fcn,varargin{:});

if exist(fuse.funconn,'dir')
    [pa,fi,dum]=fileparts(confilename);
    confilename=fullfile(pa,fi); clear pa fi    
    [isok,msg]=movefile([fuse.funconn filesep 'CComps.mat'],[confilename '.mat'],'f');
    if ~isok
        error(msg)
    else
        rmdir(fuse.funconn,'s')
    end
end


%--------------
function txt=helptxt(fcn)

T = help(fcn);
T = strrep(T,fcn,'fcm_sourceconnectivity');
T = strrep(T,'fcm_sourceconnectivity(filtfile','fcm_sourceconnectivity(outfile,filtfile');
h1 = strfind(T,'fcm_sourceconnectivity'); h1=h1(1)-3;
he = strfind(T,'Input:')-1;
start = he+8;
stop  = strfind(T,'Output:')-1;
ref   = strfind(T,'Reference');
%good  = [h1:he start:stop ref:length(T)];

txt=['  FCM_SOURCECONNECTIVITY calculates measures of functional connectivity in source space.\n\n' ...
    T(h1:he) ...
    'Input for current configuration (' fcn '):\n' ...
    '   outfile      name of new file with source connectivity data to be created.\n' ...
    T([start:stop ref:length(T)])];

