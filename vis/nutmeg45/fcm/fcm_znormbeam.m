function beam=fcm_znormbeam(beam)
% FCM_ZNORMBEAM calculates z-score at each voxel.
% beam=fcm_znormbeam(beam)


for j=1:length(beam.s)
    P = beam.s{j};
    all0 = (P==0);
    P(all0) = NaN;
    M = nanmean(P,1);
    P = P-repmat(M,[size(P,1) 1 1]);
    S = nanstd(P,[],1);
    %S(M==0 & S==0)=1;   %avoid NaN's
    beam.s{j}=P./repmat(S,[size(P,1) 1 1]);
    beam.s{j}(all0)=0;  %0 voxels usually have no true signal, therefore we avoid artifial low z scores here. 
end

if isfield(beam,'rois')
    for j=1:length(beam.rois)
        P = beam.rois{j};
        all0 = (P==0);
        P(all0) = NaN;
        M = nanmean(P,1);
        P = P-repmat(M,[size(P,1) 1 1]);
        S = nanstd(P,[],1);
        %S(M==0 & S==0)=1;   %avoid NaN's
        beam.rois{j}=P./repmat(S,[size(P,1) 1 1]);
        beam.rois{j}(all0)=0;
    end
end

if isfield(beam,'connectionfile')
    beam=rmfield(beam,'connectionfile');
end