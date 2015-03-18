function varargout=nut_nan0nan(funhandle,varargin);
% nut_nan0nan(command,arg1,arg2,...)
% temporarily replaces NaN with 0, then performs specified operation
% gets around problem when Matlab refuses to do something properly because
% of NaNs
% TODO: optionally replace corresponding output entries with NaN


nfunargout=nargout(funhandle);
if(nfunargout<0)
    nfunargout=1;
end
varargout=cell(1,nfunargout);

for ii=1:(nargin-1)
    switch class(varargin{ii})
        case {'double','single'}
            select=find(isnan(varargin{ii}));
            if(~isempty(select))
                varargin{ii}(select)=0;
            end
    end
end

[varargout{:}]=feval(funhandle,varargin{:});