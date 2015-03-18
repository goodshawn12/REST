function CC = fcm_correctspatialleakage(CC,W)
% Tries to correct for spatial leakage by subtracting the spatial filter
% correlation from the real part of coherence or from the signal amplitude
% correlation.
%
% CC = fcm_correctspatialleakage(CC,W)

if ndims(W)>2
    if size(W,2)==1
        W=permute(W,[1 3 2]);
    else
        error('The spatial filter must be scalar.')
    end
end

R = fmcorr(W); 
if strcmpi(CC.method,'ampcorr')   % it seems that we have to square this for amplitude correlation, not sure why
    R = R.^2;
end
    
dims=size(CC.coh);
if dims(1)==size(CC.frq,1)     % Legacy compatibility
    CC.coh=permute(CC.coh,[2 3 1]); 
    dims=size(CC.coh);
end    
dims(end+1:4)=1;

if dims(1)~=dims(2)         % usually true
    
    RR = zeros(dims(1),1);
    for k1=1:dims(1)
        RR(k1) = R(CC.comps(k1,1),CC.comps(k1,2));
    end
    
    if ~strcmpi(CC.method,'ampcorr') && any( (sign(real(CC.coh(:,1,1)))~=sign(RR)) & (abs(real(CC.coh(:,1,1))>.3)) & (abs(RR)>.3) )
        % once in a while the sign seems to be flipped, which may be due to
        % the ambivalence of the dipole orientation. Try to correct here.
        disp burp
        S = [sign(mean(real(CC.coh(:,1,:)),3)) sign(RR)];
        S = prod(S,2);
        RR = RR.*S; clear S
    end
        
    for k2=1:dims(2)
        for k3=1:dims(3)
            CC.coh(:,k2,k3) = CC.coh(:,k2,k3) - RR;
        end
    end
        
else        % probably does not happen anymore
       
    for k3=1:dims(3)
        for k4=1:dims(4)
            CC.coh(:,:,k3,k4) = CC.coh(:,:,k3,k4) - R;
        end
    end
    
end
if any(abs(real(CC.coh(:)))>1), warning('fcm:WrongSpatLeakCorrection','The correction for spatial leakage may not have worked correctly. Check the result.'), end


function R = fmcorr(W)
% R = fmcorr(M)

R = W.'*W;
R = R./sqrt(diag(R)*diag(R).');