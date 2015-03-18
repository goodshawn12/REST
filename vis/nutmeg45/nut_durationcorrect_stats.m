function pcorr = nut_durationcorrect_stats(T,puncorr,bin,alpha)
% NUT_DURATIONCORRECT_STATS thresholds p-values such that they must remain significant over
% at least BIN time windows.
%
% p_corr = nut_durationcorrect_stats(statval, p_uncorr [,bin [,alpha]])
% beam   = nut_durationcorrect_stats(beam [,bin [,alpha]])

if isstruct(T)
    beam=T; 
    if nargin<3, alpha=.05; else alpha=bin; end
    if nargin<2 || isempty(puncorr), bin=3; else bin=puncorr; end
    if isfield(beam,'corr') && ~isfield(beam.corr,'p_duration_corr')
        if isfield(beam.corr,'durationthreshold'), bin=beam.corr.durationthreshold;
        else beam.corr.durationthreshold=bin;
        end     
        if isfield(beam.corr,'roi_p_uncorr')
            beam.corr.roi_p_duration_corr = nut_durationcorrect_stats(beam.rois{1},beam.corr.roi_p_uncorr,bin,alpha);
            beam.corr.p_duration_corr  = ones(size(beam.s{1}));
            for k3=1:size(beam.corr.roi_p_uncorr,3)
                beam.corr.p_duration_corr(:,:,k3) = beam.R.roi2voxel_tfm * beam.corr.roi_p_duration_corr(:,:,k3);
            end
        elseif isfield(beam.corr,'p_uncorr')
            beam.corr.p_duration_corr = nut_durationcorrect_stats(beam.s{1},beam.corr.p_uncorr,bin,alpha);
        end
    end
    if isfield(beam,'snpm') && ~isfield(beam.snpm,'p_duration_corr_pos') && ~isfield(beam.snpm,'p_duration_corr_neg')
        if isfield(beam.snpm,'durationthreshold'), bin=beam.snpm.durationthreshold;
        else beam.snpm.durationthreshold=bin;
        end            
        if isfield(beam.snpm,'roi_p_uncorr_pos')
            beam.snpm.roi_p_duration_corr_pos = nut_durationcorrect_stats(beam.snpm.roi_T,beam.snpm.roi_p_uncorr_pos,bin,alpha);
            beam.snpm.p_duration_corr_pos  = ones(size(beam.s{1}));
            for k3=1:size(beam.snpm.roi_p_uncorr_pos,3)
                beam.snpm.p_duration_corr_pos(:,:,k3) = beam.R.roi2voxel_tfm * beam.snpm.roi_p_duration_corr_pos(:,:,k3);
            end
        elseif isfield(beam.snpm,'p_uncorr_pos')
            beam.snpm.p_duration_corr_pos = nut_durationcorrect_stats(beam.snpm.T,beam.snpm.p_uncorr_pos,bin,alpha);
        end
        if isfield(beam.snpm,'roi_p_uncorr_neg')
            beam.snpm.roi_p_duration_corr_neg = nut_durationcorrect_stats(beam.snpm.roi_T,beam.snpm.roi_p_uncorr_neg,bin,alpha);
            beam.snpm.p_duration_corr_neg  = ones(size(beam.s{1}));
            for k3=1:size(beam.snpm.roi_p_uncorr_neg,3)
                beam.snpm.p_duration_corr_neg(:,:,k3) = beam.R.roi2voxel_tfm * beam.snpm.roi_p_duration_corr_neg(:,:,k3);
            end
        elseif isfield(beam.snpm,'p_uncorr_neg')
            beam.snpm.p_duration_corr_neg = nut_durationcorrect_stats(beam.snpm.T,beam.snpm.p_uncorr_neg,bin,alpha);
        end        
    end
    if isfield(beam,'ttest') && ~isfield(beam.ttest,'p_duration_corr')  
        if isfield(beam.ttest,'durationthreshold'), bin=beam.ttest.durationthreshold; 
        else beam.ttest.durationthreshold=bin;
        end
        if isfield(beam.ttest,'roi_p_uncorr')
            beam.ttest.roi_p_duration_corr = nut_durationcorrect_stats(beam.ttest.T,beam.ttest.roi_p_uncorr,bin,alpha);
            beam.ttest.p_duration_corr  = ones(size(beam.s{1}));
            for k3=1:size(beam.ttest.roi_p_uncorr,3)
                beam.ttest.p_duration_corr(:,:,k3) = beam.R.roi2voxel_tfm * beam.ttest.roi_p_duration_corr(:,:,k3);
            end
        elseif isfield(beam.ttest,'p_uncorr')
            beam.ttest.p_duration_corr = nut_durationcorrect_stats(beam.ttest.T,beam.ttest.p_uncorr,bin,alpha);
        end
    end
    pcorr = beam;
    return
end
        
if nargin<4, alpha=.05; end
if nargin<3 || isempty(bin), bin=3; end

[nv,nt,nf]=size(puncorr);
if nt<bin, error('Not enough time windows for duration correction.'), end

G = false(nv,nt,nf,bin*2);
b = 1/bin .* ones(1,bin);

A = ( (puncorr<alpha) & (T>0) );  % significant positive, non-corrected
F = filter(b,1,A,[],2);  % moving average, last timewindow of significant sequences is 1.
for k=1:bin
    G(:,:,:,k) = (F==1);                    % last timewindow of significant sequence
    F = cat(2,F(:,2:end,:),zeros(nv,1,nf)); % shift to left to get remaining timewindows of sequence
end

A = ( (puncorr<alpha) & (T<0) );  % significant negative, non-corrected
F = filter(b,1,A,[],2);  % moving average, last timewindow of significant sequences is 1.
for k=1:bin
    G(:,:,:,bin+k) = (F==1);                % last timewindow of significant sequence
    F = cat(2,F(:,2:end,:),zeros(nv,1,nf)); % shift to left to get remaining timewindows of sequence
end

G = any(G,4);

pcorr = puncorr;
pcorr(~G) = 1;
