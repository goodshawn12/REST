function nuts = nut_downsample(nuts,R)
% NUT_DOWNSAMPLE  reduces sampling rate of data in nuts structure.
%                 based on DECIMATE of MATLAB's signal toolbox, but works
%                 for multiepoch signals.
%
%   nuts = nut_downsample(nuts,reduction_factor)


datadim = size(nuts.meg.data);
nout = ceil(datadim(1)/R);
nbeg = R - (R*nout - datadim(1));

[b,a] = getMatlabsDecimateChebyshefFilterDesign(R);

DD = zeros(nout,datadim(2),datadim(3));
for e=1:datadim(3)
    for c=1:datadim(2)
        odata = filtfilt(b,a,nuts.meg.data(:,c,e));
        DD(:,c,e) = odata(nbeg:R:datadim(1));
    end
end
clear odata a b
nuts.meg.data = [];
nuts.meg.data = DD;
nuts.meg.srate = nuts.meg.srate/R;
nuts.meg.latency = nuts.meg.latency(nbeg:R:datadim(1));


%--------------
function [b,a] = getMatlabsDecimateChebyshefFilterDesign(r)
% Adapted from Matlab's DECIMATE

nfilt = 8;
rip = .05;	% passband ripple in dB
[b,a] = cheby1(nfilt, rip, .8/r);
while (abs(filtmag_db(b,a,.8/r)+rip)>1e-6)
    nfilt = nfilt - 1;
    if nfilt == 0
        break
    end
    [b,a] = cheby1(nfilt, rip, .8/r);
end
if nfilt == 0
    error('Bad Chebyshev design, likely reduction_factor is too big; try mult. decimation (R=R1*R2).')
end

%--------------------------------------------------------------------------
function H = filtmag_db(b,a,f)
%FILTMAG_DB Find filter's magnitude response in decibels at given frequency.
% Adapted from Matlab's DECIMATE

nb = length(b);
na = length(a);
top = exp(-j*(0:nb-1)*pi*f)*b(:);
bot = exp(-j*(0:na-1)*pi*f)*a(:);

H = 20*log10(abs(top/bot));