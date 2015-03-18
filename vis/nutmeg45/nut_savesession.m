function nut_savesession(savepath)
% NUT_SAVESESSION   saves nutmeg session to matlab file
% Usages:
%   nut_savesession
%   nut_savesession(savepath)

global nuts;
if nargin<1
    [filename, pathname] = uiputfile('*.mat','Save File As...');
    if isequal(filename,0) | isequal(pathname,0)
        return;
    else
        savepath = fullfile(pathname,filename);
    end
else
    [pathname,filename]=fileparts(savepath);
end

nuts.sessionfile=filename;

stemp=whos('nuts');
if stemp.bytes > 2e9
    save(savepath,'-v7.3','-struct','nuts');
else
    save(savepath,'-struct','nuts');
end

