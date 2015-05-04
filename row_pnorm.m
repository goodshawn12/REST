% take the p-norm (default: 2) for a 
% matrix row by row or column by column (default: row)

function out = row_pnorm(mat,p,orientation)

if ~exist('orientation','var') || isempty(orientation)
    orientation = 'row'; end

if ~exist('p','var') || isempty(p)
    p=2; end

switch orientation
    case 'row'
        out = (sum(abs(mat).^p,2)).^(1/p);
    case 'column'
        out = (sum(abs(mat).^p)).^(1/p);
end