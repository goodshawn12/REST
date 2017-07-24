function ent = scalp_entropy( varargin )
% scalp_entropy calculates the entropy of the scalp map for each IC in the
% EEGLab data structure given as input.
%
% Input:
% eeg: EEGLab data structure with additional field
% virtual_topography. See the readme file for the
% toolbox for details on this field.
%
% Output:
% ent: Vector containing the entropy of the virtual scalp map for each IC.
% See the readme file for the toolbox for details on the virtual scalp map.

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


arg_define(1,varargin, ...
    arg('eeg', 'EEG'));

topog = eeg.virtual_topography;

n_ics=size(eeg.icawinv,2);
ent = NaN(1, n_ics);
for i= 1:n_ics
    [y,xi]=ksdensity(topog(i,:),'width',0.1);
    d1=diff(xi); % binwidth
    y=y*d1(1)+eps;
    y=y/sum(y); % normalize so that the vector y sums to one, i.e. fulfills requirement of probability distribution
    ent(i)=-sum(y.*log(y));
end
end

