function [coh,icoh]=nut_coherence(data,frequse)
% function [coh,icoh]=nut_coherence(data,frequse)
%
% data size voxels x time x trials
% data has already gone through fft e.g.
% data = fft(hanning(fs*duration).*detrend(data_time)); 
% before input here; fs (sampling rate Hz); duration (in s)
%
% frequse is indices of frequency window to use
% e.g.    linf=linspace(0,600,size(data,2));
%         frequse=[dsearchn(linf',lofreq):dsearchn(linf',hifreq)];
%
% thanks to Adrian Guggisberg for basics of code
% 20 April, 2010


for ll=1:size(data,1)
    for mm=1:size(data,1)
        Pxx = zeros(length(frequse),size(data,3));
        Pyy = zeros(length(frequse),size(data,3));
        Pxy = zeros(length(frequse),size(data,3));
        for kk=1:size(data,3)
            IC_Xx=data(ll,:,kk)';
            IC_Yy=data(mm,:,kk)';
            IC_Xx2 = abs(IC_Xx).^2;
            IC_Yy2 = abs(IC_Yy).^2;
            IC_Xy2 = IC_Yy.*conj(IC_Xx);
            Pxx(:,kk) = Pxx(:,kk) + IC_Xx2(frequse);
            Pyy(:,kk) = Pyy(:,kk) + IC_Yy2(frequse);
            Pxy(:,kk) = Pxy(:,kk) + IC_Xy2(frequse);
        end
        Pxxs = sum(Pxx,1);
        Pyys = sum(Pyy,1);
        Pxys = sum(Pxy,1);
        ICoh_timecourse = imag(Pxys./sqrt(Pxxs.*Pyys));
        icoh(ll,mm) = sum(abs(ICoh_timecourse))/size(data,3);
        Coh_timecourse = abs(Pxys./sqrt(Pxxs.*Pyys));
        coh(ll,mm) = sum(Coh_timecourse)/size(data,3);
    end
end
