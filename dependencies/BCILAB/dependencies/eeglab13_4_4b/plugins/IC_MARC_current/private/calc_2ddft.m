function [feature_2ddfts topography bandpows midims] = calc_2ddft(varargin)
% calc_2ddft calculates the 2D Discrete Fourier Transform of the topographies 
% (transposed inverse ICA weights) for each IC in the input as described in 
% Automatic Classification of Artifactual ICA-Components for Artifact 
% Removal in EEG Signals by Irene Winkler, Stefan Haufe and Michael 
% Tangermann http://www.behavioralandbrainfunctions.com/content/7/1/30
%
% Input:
% eeg: EEGLab data structure with additional fields virtual_chanlocs and
% virtual_topography. See the readme file for the toolbox for details on
% these fields.
%
% spatially_normalized_eeg (optional): EEGLab data structure in which the 
% inverse ICA weights are standardized (the columns of icawinv have mean 
% zero and variance one). If not supplied, it is calculated. 
%
% topography (optional): Matrix with rows containing the scalp map
% (electrode activation) of each component. Each row should have mean zero
% and variance one. If not supplied, it is calculated.
%
% Output:
% feature_2ddfts: Column vector containing the logarithms of the average 
% bandpower in the highest spatial frequencies in the scalp map for each IC.
%
% topography: A matrix of standardized topographies for the ICs given in
% the input. Rows contain the topographies for the ICs.
%
% bandpows: Cell array containing the 2DDFTs for each IC topography.
%
% midims: Cell array containing the central part of the 2DDFT for each IC.
% These are the matrices from which the log average bandpowers in 
% feature_2ddfts were calculated.

% This function was built by Laura Froelich (laura.frolich@gmail.com)
% using parts of topoplot.m in EEGLab. The parts taken from topoplot.m
% create a pixelated square image of the scalp map.
% Authors of topoplot.m:  Colin Humphries & Scott Makeig
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

arg_define(1,varargin, ...
    arg('eeg',[],[],'EEGLab data structure'),...
    arg('spatially_normalized_eeg', [], [],...
    'EEG struct with ICA spatial map normalized.'),...
    arg('topography', [], [],...
    'Topography.'));

if isempty(spatially_normalized_eeg)
    eeg = spatially_normalize_icdecomp(eeg);
else
    eeg = spatially_normalized_eeg;
end

if(isempty(topography))
    topography = eeg.icawinv';
end

topography = eeg.virtual_topography;
eeg.chanlocs = eeg.virtual_chanlocs;

feature_2ddfts = zeros(1, size(eeg.icawinv, 2));
bandpows = cell(length(feature_2ddfts),1);
midims = cell(length(feature_2ddfts),1);
for i = 1:length(feature_2ddfts)
    % make pixelated map
    
    Th = cell2mat({eeg.chanlocs.theta});
    Th = pi/180*Th;
    Rd = cell2mat({eeg.chanlocs.radius});
    [x y] = pol2cart(Th,Rd);
    intrad = min(1.0,max(Rd)*1.02);
    intchans = find(x <= intrad & y <= intrad);
    intx  = x(intchans);
    inty  = y(intchans);
  plotrad = min(1.0,max(Rd)*1.02);            % default: just outside the outermost electrode location
  plotrad = max(plotrad,0.5); 
    rmax = 0.5;
    squeezefac = rmax/plotrad;
    intx = intx*squeezefac;
    inty = inty*squeezefac;
    xmin = min(-rmax,min(intx)); xmax = max(rmax,max(intx));
    ymin = min(-rmax,min(inty)); ymax = max(rmax,max(inty));
    
    GRID_SCALE = 64;
    xi = linspace(xmin,xmax,GRID_SCALE);   % x-axis description (row vector)
    yi = linspace(ymin,ymax,GRID_SCALE);   % y-axis description (row vector)
    
    [Xi,Yi,Zi] = griddata(inty,intx,topography(i,intchans),yi',xi,'v4'); % interpolate data
    % pixelated map is now in Zi
    dim = length(topography);
    
    ddft = fft2(Zi); % calculate 2d fourier transform
    bandpow = abs(ddft); % calculate the power
    
    head_diameter = 18; % in cm
    nyq = (1/2)*(8/head_diameter); % nyquist rate for 8*8 electrode array on head of diameter head_diameter
    last_freq =nyq;
    halffreqs = linspace(0, (1/2)*sqrt(dim)/head_diameter, GRID_SCALE/2+1);
    freqs = [halffreqs, -halffreqs((end-1):-1:2)];
    
    % find the columns and rows in the bandpower matrix that 
    % correspond to the desired first and last frequencies.
    last_k = find(freqs >= last_freq, 1,'first')-1;
    first_freq =2/3*last_freq;
    first_k = find(freqs > first_freq, 1,'first');
    
    % we use same frequency range for both directions of image (same in x and
    % y directions), so we only need one first_freq and one first_n
    center_image_box = bandpow([first_k:last_k (GRID_SCALE-last_k+1):(GRID_SCALE-first_k+1)], first_k:last_k);
    
    log_high_freqs=log(center_image_box); % calculate log of bandpower and average of the frequencies
    % between first_freq and last_freq
    feature_2ddft = mean(mean(log_high_freqs));
    
    feature_2ddfts(i) = feature_2ddft ;
    bandpows{i} = bandpow;
    midims{i} = center_image_box;
end
