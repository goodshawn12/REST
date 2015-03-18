function nut_data_handle(command,varargin)
% NUT_DATA_HANDLE(COMMAND,MNAME1,MNAMR2 ...)
%
%   COMMAND = {'sleep','rise','flash'}
%   MNAME# - names (string) of fields you wish to address
%
% The function is designed for handling global data in NUTMEG toolbox. It
% is presumed that the toolbox keeps one and only one global variable
% (NUTVAR), which is a structure that contains all possible variables for
% this toolbox.
% One variable for storing all the data will probably make the program
% working quicker, but the problem of memory overflow will be much more
% importaint in this case. So there's a system of 2 procedures: sending
% values to a storage ('sleep') and rising them back ('rise'). 
% The first stores a field
% in the storage (a file, much like the WORKFILE which is in use now),
% deleting it from the structure. The second loads the field back to
% the NUTVAR. That helps to optimize memory and IO handling, because
% the programmer can choice whether to keep the variable in memory or
% to send it to the storage depending on what's doing and what's expected
% to be done is each part of the program.
%
% You can also use command 'flash' to clean the storage.

global nutvar;
if(~isfield(nutvar,'workfile'))
    error(' The WORKFILE name isn''t specified');
end

switch(lower(command))
    case{'flash'}
        try
            delete(nutvar.workfile);
        catch 
            disp(['Creating file ' nutvar.workfile]);
        end
        a = []; % A silly way to create an empty .mat file.
        save(nutvar.workfile,'a');
    case{'rise'}
        for(Q=1:(nargin-1))
            if(isfield(nutvar,varargin{Q})); continue; end;
            try
                temp = load(nutvar.workfile,varargin{Q});
                nutvar = setfield(nutvar,varargin{Q},getfield(temp,varargin{Q}));
            catch
                error('NUT_DATA_HANDLE Warning: You have addressed (rise) an unknown global variable');
            end
        end
    case{'sleep'}    
        for(Q=1:(nargin-1))
            try
                eval([varargin{Q} '=nutvar.' varargin{Q} ';']);
                save(nutvar.workfile,varargin{Q},'-APPEND');
                nutvar = rmfield(nutvar,varargin{Q});
            catch                
                list = who('-file',nutvar.workfile);
                if(~ismember(varargin(Q),list))
                    error('NUT_DATA_HANDLE Warning: You have addressed (sleep) an unknown global variable');
                else
                    disp('Warning: narcotizing already sleeping variable. Don''t waste machine time!');
                end
            end
        end
    otherwise
        error('Unknown command');
end
return;