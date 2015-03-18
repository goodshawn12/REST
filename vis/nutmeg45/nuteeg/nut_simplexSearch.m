%%
% Simplex search algorithm based on Nelder & Mead, 1965.
% @param funfcn the cost function to minimize.
% @x the initiallization point. Should contain 3 dimensions. 1st two
% dimensions describe the search point where convergence is defined when
% the sum of the distances between corresponding columns is less than tolX.
% 3rd dimension contains the different initiallization points.
% @options parameters for the search.
% @author Daniel D.E. Wong
%%
function [x,y,cflag] = nut_simplexSearch(funfcn,x,options,varargin)

if size(x,3) ~= size(x,1)*size(x,2)+1
    error('Invalid dimensions. n+1 simplex points required.');
    return;
end

if ~exist('options') | ~isfield(options,'tolX'); options.tolX = 1; end;
if ~exist('options') | ~isfield(options,'maxIt'); options.maxIt = 1000; end;
if ~exist('options') | ~isfield(options,'alpha'); options.alpha = 1; end;
if ~exist('options') | ~isfield(options,'gamma'); options.gamma = 2; end;
if ~exist('options') | ~isfield(options,'rho'); options.rho = 0.5; end;
if ~exist('options') | ~isfield(options,'sigma'); options.sigma = 0.5; end;

y = zeros(size(x,3),1);
for i = 1:size(x,3)
    y(i) = funfcn(x(:,:,i),varargin{:});
end

cflag = 0;
it = -1;
while 1
    it = it + 1;
    if it > options.maxIt; warning('Maximum iterations exceeded'); break; end;
    
    % Sort
    [y,yidx] = sort(y);
    x = x(:,:,yidx);
    
    % Check convergence
    x0 = mean(x(:,:,1:end),3);
    cflag = 1;
    for i = 1:size(x,3)
        for j = 1:size(x,1)
            if norm(x(j,:,i)-x0(j,:)) > options.tolX
                cflag = 0;
                break;
            end
        end
        if ~cflag; break; end;
    end
    if cflag; break; end;
    
    x0 = mean(x(:,:,1:end-1),3);
    
    % Reflection
    xr = x0 + options.alpha*(x0-x(:,:,end));
    yr = funfcn(xr,varargin{:});
    if y(1) <= yr & yr < y(end-1)
        x(:,:,end) = xr;
        y(end) = yr;
        continue;
    end
    
    % Expansion
    if yr < y(1)
        xe = x0 + options.gamma*(x0-x(:,:,end));
        ye = funfcn(xe,varargin{:});
        if ye < yr
            x(:,:,end) = xe;
            y(end) = ye;
        else
            x(:,:,end) = xr;
            y(end) = yr;
        end
        continue;
    end
    
    % Contraction
    if yr >= y(end-1)
        xc = x(:,:,end) + options.rho*(x0-x(:,:,end));
        yc = funfcn(xc,varargin{:});
        if yc < y(end)
            x(:,:,end) = xc;
            y(end) = yc;
            continue;
        end
    end
    
    % Shrink
    x1 = repmat(x(:,:,1),[1 1 size(x,3)]);
    x = x1 + options.sigma*(x-x1);
    continue;
end
x = x0;
y = funfcn(x0,varargin{:});