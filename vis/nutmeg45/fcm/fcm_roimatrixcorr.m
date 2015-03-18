function MM=fcm_roimatrixcorr(MM,clindata)
% MM=fcm_roimatrixcorr(MM,clindata)

FDR = 0.1;

[numpat,numroi,dum,nt,nf]=size(MM.coh);
if length(clindata)~=numpat,  error('Behavioral data and connectivity data do not have the same number of subjects.'), end

clindata = clindata(:);

MM.r = zeros(numroi,numroi,nt,nf);
MM.p_uncorr = ones(numroi,numroi,nt,nf);
for rr = 1:numroi
    for ff = 1:nf
        for tt = 1:nt
            [MM.r(:,rr,tt,ff),MM.p_uncorr(:,rr,tt,ff)] = corr(MM.coh(:,:,rr,tt,ff),clindata,'rows','pairwise');
        end
    end
end

uti = isfinite(tril(nan(numroi),0));

MM.p_FDR_corr = ones(size(MM.p_uncorr));
MM.cutoff = zeros(nt,nf);
for ff=1:nf
    for tt=1:nt
        temp=nut_FDR(MM.p_uncorr(:,:,nt,nf),FDR);
        [out,MM.cutoff(tt,ff)] = nut_FDR(temp(uti),FDR);
        temp(:,:)=0;
        temp(uti)=out;
        temp = temp + temp.' + diag(ones(1,numroi));
        MM.p_FDR_corr(:,:,tt,ff) = temp;
    end
end
MM.info='Correlation';