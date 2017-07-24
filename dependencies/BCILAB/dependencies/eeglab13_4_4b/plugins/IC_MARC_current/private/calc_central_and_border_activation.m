function [f_cent,f_border] = calc_central_and_border_activation(X,chanlocs,transform)
% Calculate the Central Activation and Border activation features described
% in Automatic Classification of Artifactual ICA-Components for Artifact 
% Removal in EEG Signals by Irene Winkler, Stefan Haufe, and Michael
% Tangermann
% http://www.behavioralandbrainfunctions.com/content/7/1/30
%
% In:
%   X : scalp map vector (elsewhere referred to as the topography of an IC,
%       the columns of the inverse ICA weight matrix)
%
%   Chanlocs: channel locations
%
%   Transform : optional transformation ensuring that +Z is up and +X is the nose direction (default: eye(3))
%
% Out:
%   c : central activation
%   b : border activation

% Authors: Christian Kothe and Laura Froelich (laura.frolich@gmail.com)
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

if ~exist('transform','var')
    transform = eye(3); end


% look up channel locations if necessary
if ~all(isfield(chanlocs,{'X','Y','Z'})) || all(cellfun('isempty',{chanlocs.X}))
    chanlocs = pop_chanedit(chanlocs,'lookup','standard-10-5-cap385.elp'); end

% find the admissible channel indices and map channels into seamless spherical coordinate system
ok = find(~(cellfun('isempty',{chanlocs.X}) | cellfun('isempty',{chanlocs.Y}) | cellfun('isempty',{chanlocs.Z})));
chanlocs = chanlocs(ok); X = X(ok);

% transform channel locations...
transformed = [[chanlocs.X]' [chanlocs.Y]' [chanlocs.Z]']*eye(3);
tmp = num2cell(transformed(:,1)); [chanlocs.X] = deal(tmp{:});
tmp = num2cell(transformed(:,2)); [chanlocs.Y] = deal(tmp{:});
tmp = num2cell(transformed(:,3)); [chanlocs.Z] = deal(tmp{:});

% get the cz channel index
[maxval,cz] = max([chanlocs.Z]);

% find the central electrode group
[px,py] = cart2sph([chanlocs.Z],[chanlocs.X],[chanlocs.Y]);
distances = hypot(px(cz)-px,py(cz)-py);

% constants
central_group_size = min(13, size(chanlocs,2));
outermost_group_size = min(18, size(chanlocs,2));
occipital_group_size = min(17, size(chanlocs,2));

[dummy,indices] = sort(distances,'ascend');
central_group = indices(1:central_group_size);

% compute central activation
f_cent = log(mean(abs(X(central_group))));

% find the k outermost electrodes
[dummy,indices] = sort(distances,'descend');
outer_group = indices(1:outermost_group_size);

% compute the border-activation feature
f_border = -1;
[dummy,maxidx] = max(abs(X));
if ismember(maxidx,outer_group)
    f_border = 1; % if the most active electrode in the pattern is part
    % of the outermost group, set f_border to 1
else
    % compute the 5 frontal sectors...
    angles = rad2deg(atan2(py-py(cz),px-px(cz)));
    for sector = -90:36:54
        % channel subset
        chanmask = setdiff(find(within(angles,sector,sector+36)),cz);
        if ~isempty(chanmask)
            % maximum value in local group
            [maxval,imaxval] = max(abs(X(chanmask)));
            % distances in local group
            [maxdist,imaxdist] = max(distances(chanmask));
            mu = mean(X(chanmask));
            st = std(X(chanmask));
            if imaxval == imaxdist && (maxval - mu)/st >= 2
                f_border = 1; end
        end
    end
end
% compute the occipital sector...
[dummy,oz] = min([chanlocs.X]);
ozdistances = hypot(px(oz)-px,py(oz)-py);
[dummy,indices] = sort(ozdistances,'ascend');
occipital_group = indices(1:occipital_group_size);
% maximum value in local group
[maxval,imaxval] = max(abs(X(occipital_group)));
% distances in local group
[maxdist,imaxdist] = max(distances(occipital_group));
mu = mean(X(occipital_group));
st = std(X(occipital_group));
if imaxval == imaxdist && (maxval - mu)/st >= 2
    f_border = 1; end
end

function tf = within(x,a,b)
if b<a b = b+360; end
tf = (x>=a & x<b) | ((x+360)>=a & (x+360)<b) | ((x-360)>=a & (x-360)<b);
end

function x = rad2deg(x)
x = (180/pi) * x;
end
