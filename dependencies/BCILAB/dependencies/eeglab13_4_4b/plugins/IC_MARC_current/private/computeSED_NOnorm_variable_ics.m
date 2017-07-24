function [out,medie_left,medie_right]=computeSED_NOnorm_variable_ics(topog,chanlocs,n, n_ics)
% computeSED_NOnorm() - Computes Spatial Eye Difference feature
% without normalization 
%
% Usage:
%   >> [out,medie_left,medie_right]=computeSED_NOnorm(topog,chanlocs,n);
%
% Inputs:
%   topog      - topographies vector
%   chanlocs   - EEG.chanlocs struct
%   n          - number of electrodes
%
% Outputs:
%   out        - SED values
%   medie_left - Left Eye area average values
%   medie_right- Right Eye area average values
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
%
% 23/4/2013: modified by Laura Froelich to accomodate a number of ic's not equal to the number
% of electrodes

% Define scalp zones

% Find electrodes in Left Eye area (LE)
dimleft=0; %number of LE electrodes
index1=zeros(1,n); %indexes of LE electrodes

for k=1:n 
    if (-61<chanlocs(1,k).theta) && (chanlocs(1,k).theta<-35) && (chanlocs(1,k).radius>0.30) %electrodes are in LE
        dimleft=dimleft+1; %count electrodes
        index1(1,dimleft)=k; 
    end
end
    
 % Find electrodes in Right Eye area (RE)
    dimright=0; %number of RE electrodes
    index2=zeros(1,n); %indexes of RE electrodes
    for g=1:n 
        if (34<chanlocs(1,g).theta) && (chanlocs(1,g).theta<61) && (chanlocs(1,g).radius>0.30) %electrodes are in RE
            dimright=dimright+1; %count electrodes
            index2(1,dimright)=g; 
        end
    end
    
if(dimleft*dimright==0)
        disp('ERROR: no channels included in some scalp areas.')
        disp('Check channels distribution and/or change scalp areas definitions in computeSAD.m and computeSED_NOnorm.m')
        disp('ADJUST session aborted.')
        return
end
% Outputs

     out=zeros(n_ics,1); %memorizes SED
      medie_left=zeros(n_ics,1); %memorizes LE mean value
      medie_right=zeros(n_ics,1); %memorizes RE mean value
      
% Output computation

for i=1:n_ics  % for each topography
 %create LE electrodes vector
    left=zeros(1,dimleft);
    for h=1:dimleft
        left(1,h)=topog(i,index1(1,h));
    end
    
   %create RE electrodes vector
    right=zeros(1,dimright);
    for h=1:dimright
        right(1,h)=topog(i,index2(1,h));
    end
    
    %compute features
    out(i)=abs(mean(left)-mean(right));% SED not normalized
    medie_left(i)=mean(left);
    medie_right(i)=mean(right);
    
    
end

   
