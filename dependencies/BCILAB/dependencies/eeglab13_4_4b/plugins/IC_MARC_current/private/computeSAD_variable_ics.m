function [rapp,var_front,var_back,mean_front,mean_back]=computeSAD_variable_ics(topog,chanlocs,n, n_ics)
% computeSAD() - Computes Spatial Average Difference feature 
%
% Usage:
%   >> [rapp,var_front,var_back,mean_front,mean_back]=computeSAD(topog,chanlocs,n, n_ics);
%
% Inputs:
%   topog      - matrix of topographies
%   chanlocs   - EEG.chanlocs struct
%   n          - number of electrodes
%   n_ics      - number of ICs
%
% Outputs:
%   rapp       - SAD values
%   var_front  - Frontal Area variance values
%   var_back   - Posterior Area variance values
%   mean_front - Frontal Area average values
%   mean_back  - Posterior Area average values
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
% 22/4/2013: Modified by Laura Froelich. Now accomodates a number of ICs 
% different from the number of channels

% Define scalp zones

% Find electrodes in Frontal Area (FA)
dimfront=0; %number of FA electrodes
index1=zeros(1,n); %indexes of FA electrodes
for k=1:n
    if (abs(chanlocs(1,k).theta)<60) && (chanlocs(1,k).radius>0.40) %electrodes are in FA
        dimfront=dimfront+1; %count electrodes
        index1(1,dimfront)=k; 
    end
end

 % Find electrodes in Posterior Area (PA)
    dimback=0;
    index3=zeros(1,n);
    for h=1:n 
        if (abs(chanlocs(1,h).theta)>110) 
            dimback=dimback+1; 
            index3(1,dimback)=h; 
        end
    end
 
    if dimfront*dimback==0
        disp('ERROR: no channels included in some scalp areas.')
        disp('Check channels distribution and/or change scalp areas definitions in computeSAD.m and computeSED_NOnorm.m')
        disp('ADJUST session aborted.')
        return
    end
    
% Outputs

     rapp=zeros(1,n_ics); % SAD
      mean_front=zeros(n_ics,1); % FA electrodes mean value
      mean_back=zeros(n_ics,1); % PA electrodes mean value
      var_front=zeros(n_ics,1); % FA electrodes variance value
      var_back=zeros(n_ics,1); % PA electrodes variance value

% Output computation

for i=1:n_ics % for each topography
    
 %create FA electrodes vector
    front=zeros(1,dimfront);
    for h=1:dimfront
        front(1,h)=topog(i,index1(1,h));
    end
    
  %create PA electrodes vector
    back=zeros(1,dimback);
    for h=1:dimback
        back(1,h)=topog(i,index3(1,h));
    end
    
   %compute features
    rapp(i)=abs(mean(front))-abs(mean(back)); % SAD
    mean_front(i)=mean(front);
    mean_back(i)=mean(back);
    var_back(i)=var(back);
    var_front(i)=var(front);
end

