function icaacts = temporally_normalize_icdecomp(icaacts)
% Standardizes IC time series to have mean zero and variance one.
%
% Input:
% icaacts: IC activations.

% Output:
% icaact: Standardized IC time series 
% (rows of icaact have mean zero and variance one).

% Copyright (C) 2013  Laura Froelich (laura.frolich@gmail.com)
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

if length(size(icaacts))==3
    nepochs = size(icaacts,3);
    for iepoch=1:nepochs
        icaacts(:,:,iepoch) = zscore(icaacts(:,:,iepoch),[],2);
    end
else
icaacts = zscore(icaacts,[],2);
end
