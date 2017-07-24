function discon_meas = localized_discontinuity_measure(topographies, channel_locations, n_ics)
% localized_discontinuity_measure calculates a measure of how discontinuous
% an IC scalp map is. A high value indicates that the IC might represent an
% artifact. A very similar measure was used for the ADJUST toolbox,
% http://www.unicog.org/pm/pmwiki.php/MEG/RemovingArtifactsWithADJUST 
% Our weighting scheme for electrodes surrounding the electrode for which
% the measure is being calculated is a little different.
%
% Input:
% topographies: Matrix of topographies (scalp maps) of ICs. Each row
% corresponds to an IC, and each column to an electrode, so it is the
% transpose of the field icawinv in EEGLab data structures.
%
% channel_locations: Channel locations of the electrodes whose activations
% are given in topographies.
%
% n_ics: The number of ICs represented in topographies.
%
% Output:
% discon_meas: A vector of discontinuity measures for each IC. For each IC,
% the discontinuity measure is found by first calculating how peaked the
% spatial map is for each electrode. The discontinuity measure is then the
% highest value found.

% This code is based heavily on the code from the ADJUST toolbox. It was
% modified by Laura Froelich 24/4/2013.
%
% Copyright (C) 2009 Andrea Mognon and Marco Buiatti, 
% Center for Mind/Brain Sciences, University of Trento, Italy
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

xpos=[channel_locations.X];ypos=[channel_locations.Y];zpos=[channel_locations.Z];
pos=[xpos',ypos',zpos'];
discon_meas=zeros(1,n_ics);

for ic=1:n_ics
    aux=zeros(1,length(channel_locations));
    for el=1:length(channel_locations)-1
        
        P=pos(el,:); %position of current electrode
        d=pos-repmat(P,length(channel_locations),1);
        dist=sqrt(sum((d.*d),2));
        
        [y,I]=sort(dist, 'ascend'); % sort electrodes in ascending order of distance, so closest are first in list
        repchas=I(2:end);
        weightdist=exp(-y(2:end)); % weights computed wrt distance
        weightdist=weightdist/sum(weightdist); % ensure that weights sum to one
        
        aux(el)=abs(topographies(ic,el)-sum(weightdist.*topographies(ic,repchas)')); % Since the sum
        % of weights is one, the quantity
        % sum(weightdist.*topographies(ic,repchas)'), is a weighted
        % average. Hence this calculation finds the difference between the
        % activation in the current electrode and the average activation
        % surrounding it.
    end
    
    discon_meas(ic)=max(aux);
end

