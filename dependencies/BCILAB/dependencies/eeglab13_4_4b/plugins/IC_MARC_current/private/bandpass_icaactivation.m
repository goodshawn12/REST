function eeg = bandpass_icaactivation(eeg, rp, rs, locutoff, hicutoff)
% bandpass_icaactivation applies a bandpass filter to each IC time series 
% in the field icaact_filtered_resampled,
% to each epoch separately.
% This function is built from part of iirfilt.m in EEGLab.
%
% Input:
% eeg: EEGLab data structure. The data structure must have an ica
% decomposition as well as the additional fields icaact_filtered_resampled
% and virtual_srate. See the readme file for the toolbox for details on
% these fields.
%
% rp: ripple amplitude in dB in the pass band (rp: ripple pass)
%
% rs: ripple amplitude in dB in the stop band (rs: ripple stop)
%
% locutoff: lower limit of the pass band
%
% hicutoff: upper limit of the pass band
%
% Output:
% eeg: The EEGLab data structure given as input, but with filtered time
% series in the field icaact_filtered_resampled.

% This function was built by Laura Froelich (laura.frolich@gmail.com)
% from part of iirfilt.m in EEGLab.
% Authors of iirfilt.m: Maksym Pozdin (mpozdin.ece04@gtalumni.org, IOL/ONRC,2004), 
%          with Arnaud Delorme and Scott Makeig (SCCN/INC/UCSD, La Jolla CA)
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

icaact=eeg.icaact_filtered_resampled;
icaact = double(icaact);
nepochs = size(icaact,3);
n_ics = size(icaact,1);
    
trans_bw = 1;
trans_bw_low = 1/2;
nyq = eeg.virtual_srate/2;

% highpass filter
    ws=(locutoff-trans_bw_low)/nyq;
    wp=(locutoff)/nyq;
    [N,wn] = ellipord(wp,ws,rp,rs);
    
    [b,a]=ellip(N,rp,rs,wn,'high');
    
    for ic=1:n_ics
        for epoch=1:nepochs
            icaact(ic,:, epoch) = filtfilt(b,a,icaact(ic,:, epoch));
        end
    end
    
%  lowpass filter

    ws=(hicutoff+trans_bw)/nyq;
    wp=(hicutoff)/nyq;
    [N,wn] = ellipord(wp,ws,rp,rs);
    
    [b,a]=ellip(N,rp,rs,wn);

    for ic=1:n_ics
        for epoch=1:nepochs
            icaact(ic,:, epoch) = filtfilt(b,a,icaact(ic,:, epoch));
        end
    end
eeg.icaact_filtered_resampled = icaact;
end

