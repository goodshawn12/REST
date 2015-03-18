function SB=fcm_meansensorcomcoh(SC)

CFT=fcm_fisherZtf(SC.coh,'ccohere');

sl=unique(SC.comps(:))';
numsens=length(sl);
num2   =size(SC.coh,2);
temp = zeros(numsens,num2);

for k=sl
    f = any(SC.comps==k,2);
    temp(k,:) = mean(CFT(f,:),1);
end

% Inverse Fisher
temp = fcm_fisherZtf(temp,'ccohere',-1);

% Output structure
SB=struct('s',{{[]}},'frq',SC.frq);    

SB.s{1} = imag( temp );
SB.s{2} = abs(temp).^2;
SB.s{3} = real( temp ); 