function [p, h, stats] = signrank(x,y,varargin)
%SIGNRANK Wilcoxon signed rank test for zero median, fixed by nutmeg to actually give you the damned sign
%   P = SIGNRANK(X) performs a two-sided signed rank test of the hypothesis
%   that the data in the vector X come from a distribution whose median
%   (and mean, if it exists) is zero, and returns the p-value from the
%   test.  P is the probability of observing the given result, or one more
%   extreme, by chance if the null hypothesis ("median is zero") is true.
%   Small values of P cast doubt on the validity of the null hypothesis.
%   The data are assumed to come from a continuous distribution, symmetric
%   about its median.
%
%   P = SIGNRANK(X,M) performs a two-sided test of the hypothesis that the
%   data in the vector X come from a distribution whose median is M.  M
%   must be a scalar.
%
%   P = SIGNRANK(X,Y) performs a paired, two-sided test of the hypothesis
%   that the difference between the matched samples in the vectors X and Y
%   comes from a distribution whose median is zero.  The differences X-Y
%   are assumed to come from a continuous distribution, symmetric about its
%   median.  X and Y must be the same length.  The two-sided p-value is
%   computed by doubling the most significant one-sided value.
%
%   SIGNRANK treats NaNs in X or Y as missing values, and removes them.
%
%   [P,H] = SIGNRANK(...) returns the result of the hypothesis test,
%   performed at the 0.05 significance level, in H.  H=0 indicates that
%   the null hypothesis ("median is zero") cannot be rejected at the 5%
%   level. H=1 indicates that the null hypothesis can be rejected at the
%   5% level.
%
%   [P,H] = SIGNRANK(...,'alpha',ALPHA) returns the result of the hypothesis
%   test performed at the significance level ALPHA.
%
%   [P,H] = SIGNRANK(...,'method',METHOD) computes the p-value using an
%   exact algorithm if METHOD is 'exact', or a normal approximation if
%   METHOD is 'approximate'.  The default is to use an exact method for
%   small samples.
%
%   [P,H,STATS] = SIGNRANK(...) returns STATS, a structure with one or two
%   fields.  The field 'signedrank' contains the value of the signed rank
%   statistic.  If P is calculated using a normal approximation, then the
%   field 'zval' contains the value of the normal (Z) statistic.
%
%   See also SIGNTEST, RANKSUM, TTEST, ZTEST.

%   References:
%      [1] Hollander, M. and D. A. Wolfe.  Nonparametric Statistical
%          Methods. Wiley, 1973.
%      [2] Gibbons, J.D.  Nonparametric Statistical Inference,
%          2nd ed.  M. Dekker, 1985.

%   Copyright 1993-2005 The MathWorks, Inc. 
%   $Revision: 1.16.4.4 $  $Date: 2005/06/21 19:46:36 $
%   Modified by SSD to yield the f'in sign like even your mother
%     would tell you is fairly obligatory.

% Check most of the inputs now
alpha = 0.05;
if nargin>2 && isnumeric(varargin{1})
   % Grandfathered syntax:  signrank(x,y,alpha)
   alpha = varargin{1};
   varargin(1) = [];
end
oknames = {'alpha' 'method'};
dflts   = {alpha   ''};
[eid,emsg,alpha,method] = statgetargs(oknames,dflts,varargin{:});
if ~isempty(eid)
   error(sprintf('stats:signrank:%s',eid),emsg);
end

if ~isscalar(alpha)
   error('stats:signrank:BadAlpha','SIGNRANK requires a scalar ALPHA value.');
end
if ~isnumeric(alpha) || isnan(alpha) || (alpha <= 0) || (alpha >= 1)
   error('stats:signrank:BadAlpha','SIGNRANK requires 0 < ALPHA < 1.');
end

if nargin < 2 || isempty(y)
    y = zeros(size(x));
elseif isscalar(y)
    y = repmat(y, size(x));
end

if ~isvector(x) || ~isvector(y)
    error('stats:signrank:InvalidData',...
          'SIGNRANK requires vector rather than matrix data.');
elseif numel(x) ~= numel(y)
    error('stats:signrank:InputSizeMismatch',...
    'SIGNRANK requires the data vectors to have the same number of elements.');
end

diffxy = x(:) - y(:);

% Remove missing data
diffxy(isnan(diffxy)) = [];
if (length(diffxy)==0)
   error('stats:signrank:NotEnoughData','No data remaining after removal of NaNs.');
end

nodiff = find(diffxy == 0);
diffxy(nodiff) = [];
n = length(diffxy);

if (n == 0)         % degenerate case, all ties
    p = 1;
    if (nargout > 1)
        h = 0;
        if (nargout > 2)
            stats.signedrank = 0;
        end
    end
    return
end

% Now deal with the method argument
if isempty(method)
   if n<=15
      method = 'exact';
   else
      method = 'approximate';
   end
elseif ischar(method)
   okmethods = {'exact' 'approximate' 'oldexact'};
   j = strmatch(lower(method),okmethods);
   if isempty(j)
      error('stats:signrank:BadMethod',...
            'METHOD must be ''exact'' or ''approximate''.');
   end
   method = okmethods{j};
else
   error('stats:signrank:BadMethod',...
         'METHOD must be ''exact'' or ''approximate''.');
end

% Find negative differences and ranks of absolute differences
neg = find(diffxy<0);
[tierank, tieadj] = tiedrank(abs(diffxy));

% Compute signed rank statistic (most extreme version)
w = sum(tierank(neg));
if(w < n*(n+1)/2-w)
    sgn=-1;
else
    sgn=1;
end
w = min(w, n*(n+1)/2-w);

if isequal(method,'approximate')
    z = (w-n*(n+1)/4) / sqrt((n*(n+1)*(2*n+1) - tieadj)/24);
    p = 2*normcdf(z,0,1);
    if (nargout > 2)
        stats.zval = z*sgn;
    end
elseif isequal(method,'oldexact')
    % Enumerates all possibilities and does not adjust for ties
    allposs = (ff2n(n))';
    idx = (1:n)';
    idx = idx(:,ones(2.^n,1));
    pranks = sum(allposs.*idx,1);
    tail = 2*length(find(pranks <= w)); % two side.

    % Avoid p>1 if w is in the middle and is double-counted
    p = min(1, tail./(2.^n));
else % isequal(method,'exact')
    p = statsrexact(tierank,w);
    p = min(1, 2*p);   % two-sided, don't double-count the middle value
end

if nargout > 1
    h = (p<=alpha);
    if (nargout > 2)
        stats.signedrank = w*sgn;
    end
end





function [eid,emsg,varargout]=statgetargs(pnames,dflts,varargin)
%STATGETARGS Process parameter name/value pairs for statistics functions
%   [EID,EMSG,A,B,...]=STATGETARGS(PNAMES,DFLTS,'NAME1',VAL1,'NAME2',VAL2,...)
%   accepts a cell array PNAMES of valid parameter names, a cell array
%   DFLTS of default values for the parameters named in PNAMES, and
%   additional parameter name/value pairs.  Returns parameter values A,B,...
%   in the same order as the names in PNAMES.  Outputs corresponding to
%   entries in PNAMES that are not specified in the name/value pairs are
%   set to the corresponding value from DFLTS.  If nargout is equal to
%   length(PNAMES)+1, then unrecognized name/value pairs are an error.  If
%   nargout is equal to length(PNAMES)+2, then all unrecognized name/value
%   pairs are returned in a single cell array following any other outputs.
%
%   EID and EMSG are empty if the arguments are valid.  If an error occurs,
%   EMSG is the text of an error message and EID is the final component
%   of an error message id.  STATGETARGS does not actually throw any errors,
%   but rather returns EID and EMSG so that the caller may throw the error.
%   Outputs will be partially processed after an error occurs.
%
%   This utility is used by some Statistics Toolbox functions to process
%   name/value pair arguments.
%
%   Example:
%       pnames = {'color' 'linestyle', 'linewidth'}
%       dflts  = {    'r'         '_'          '1'}
%       varargin = {{'linew' 2 'nonesuch' [1 2 3] 'linestyle' ':'}
%       [eid,emsg,c,ls,lw] = statgetargs(pnames,dflts,varargin{:})    % error
%       [eid,emsg,c,ls,lw,ur] = statgetargs(pnames,dflts,varargin{:}) % ok

%   Copyright 1993-2004 The MathWorks, Inc. 
%   $Revision: 1.4.2.1 $  $Date: 2003/11/01 04:28:41 $ 

% We always create (nparams+2) outputs:
%    one each for emsg and eid
%    nparams varargs for values corresponding to names in pnames
% If they ask for one more (nargout == nparams+3), it's for unrecognized
% names/values

% Initialize some variables
emsg = '';
eid = '';
nparams = length(pnames);
varargout = dflts;
unrecog = {};
nargs = length(varargin);

% Must have name/value pairs
if mod(nargs,2)~=0
    eid = 'WrongNumberArgs';
    emsg = 'Wrong number of arguments.';
else
    % Process name/value pairs
    for j=1:2:nargs
        pname = varargin{j};
        if ~ischar(pname)
            eid = 'BadParamName';
            emsg = 'Parameter name must be text.';
            break;
        end
        i = strmatch(lower(pname),pnames);
        if isempty(i)
            % if they've asked to get back unrecognized names/values, add this
            % one to the list
            if nargout > nparams+2
                unrecog((end+1):(end+2)) = {varargin{j} varargin{j+1}};
                
                % otherwise, it's an error
            else
                eid = 'BadParamName';
                emsg = sprintf('Invalid parameter name:  %s.',pname);
                break;
            end
        elseif length(i)>1
            eid = 'BadParamName';
            emsg = sprintf('Ambiguous parameter name:  %s.',pname);
            break;
        else
            varargout{i} = varargin{j+1};
        end
    end
end

varargout{nparams+1} = unrecog;
