function nut_cleanup(TODO,strings)
% NUT_CLEANUP  searches folder and all subdirectories for files and deletes them.
%
%  nut_cleanup(directories,files)
%     directories    names of one or several folders which should be cleaned,
%                    including all the subdirectories of these folders!
%     files          filenames to be cleaned. This can include * as
%                    wildcards, e.g., {'cov*.mat' '*_spatnorm.mat'} 


warning('off','MATLAB:DELETE:FileNotFound')
currdir=pwd;

if ~iscell(TODO) && ~isempty(TODO), TODO={TODO}; end
if ~iscell(strings), strings={strings}; end

while ~isempty(TODO)
    newtodo = [];
    for k=1:length(TODO)
        if isdir(TODO{k})
            cd(TODO{k})
            del(strings);
            d = dir('*'); 
            d(1:2)=[]; % remove parent directories (. and ..)
            subdirs = {d([d.isdir]).name};
            subdirs = cellfun(@(x) [pwd filesep x],subdirs,'UniformOutput',false); % convert to full path
            newtodo = [newtodo subdirs];
            cd ..
        end
    end
    TODO = newtodo;
end

cd(currdir);
warning('on','MATLAB:DELETE:FileNotFound')

%----------------
function del(strings)

for k=1:length(strings)
    delete(strings{k});
end