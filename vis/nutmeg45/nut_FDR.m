function [P,u,oldout] = nut_FDR(uP,q,old)
% Calculates False Discovery Rate
%   [p_corr,cutoff] = nut_FDR(p_uncorr,q)
%   beam = nut_FDR(beam,q)
%
% p_uncorr  - Uncorrected p-values 
% q         - critical expected False Discovery Rate
%
% cutoff    - signifance threshold of uncorrected p-values.
% p_corr    - FDR corrected p-value.
% 
%___________________________________________________________________________
%
% The Benjamini & Hochberch (1995) False Discovery Rate (FDR) procedure
% finds a threshold u such that the expected FDR is at most q.  
%
% Benjamini & Hochberg (1995), "Controlling the False Discovery Rate: A
% Practical and Powerful Approach to Multiple Testing". J Royal Stat Soc,
% Ser B.  57:289-300.
%___________________________________________________________________________
% Based on snpm_uc_FDR.m and snpm_P_FDR.m from the SnPM3 toolbox:
%    http://www.fil.ion.ucl.ac.uk/spm/snpm

if nargin>2         % legacy compatibility
    if isempty(q)
        q=old; clear old
    else
        help nut_FDR
        error('This function does not support baseline windows anymore.')
    end
end

if isstruct(uP)
    beam=uP;
    if isfield(beam,'corr') && ~isfield(beam.corr,'p_FDR_corr')        
        if isfield(beam.corr,'FDR'), q=beam.corr.FDR;
        else beam.corr.FDR=q;
        end                 
        if isfield(beam.corr,'roi_p_uncorr')
            [beam.corr.roi_p_duration_corr,beam.corr.cutoff] = nut_FDR(beam.corr.roi_p_uncorr,q);
            beam.corr.p_FDR_corr  = ones(size(beam.s{1}));
            for k3=1:size(beam.corr.roi_p_uncorr,3)
                beam.corr.p_FDR_corr(:,:,k3) = beam.R.roi2voxel_tfm * beam.corr.roi_p_duration_corr(:,:,k3);
            end
        elseif isfield(beam.corr,'p_uncorr')
            [beam.corr.p_FDR_corr,beam.corr.cutoff] = nut_FDR(beam.corr.p_uncorr,q);
        end
    end
    if isfield(beam,'snpm') && ~isfield(beam.snpm,'p_FDR_corr_pos') && ~isfield(beam.snpm,'p_FDR_corr_neg')        
        if isfield(beam.snpm,'FDR'), q=beam.snpm.FDR;
        else beam.snpm.FDR=q;
        end                   
        if isfield(beam.snpm,'roi_p_uncorr_pos')
            [beam.snpm.roi_p_duration_corr_pos,beam.snpm.cutoff_pos] = nut_FDR(beam.snpm.roi_p_uncorr_pos,q);
            beam.snpm.p_FDR_corr_pos  = ones(size(beam.s{1}));
            for k3=1:size(beam.snpm.roi_p_uncorr_pos,3)
                beam.snpm.p_FDR_corr_pos(:,:,k3) = beam.R.roi2voxel_tfm * beam.snpm.roi_p_duration_corr_pos(:,:,k3);
            end
        elseif isfield(beam.snpm,'p_uncorr_pos')
            [beam.snpm.p_FDR_corr_pos,beam.snpm.cutoff_pos] = nut_FDR(beam.snpm.p_uncorr_pos,q);
        end
        if isfield(beam.snpm,'roi_p_uncorr_neg')
            [beam.snpm.roi_p_duration_corr_neg,beam.snpm.cutoff_neg] = nut_FDR(beam.snpm.roi_p_uncorr_neg,q);
            beam.snpm.p_FDR_corr_neg  = ones(size(beam.s{1}));
            for k3=1:size(beam.snpm.roi_p_uncorr_neg,3)
                beam.snpm.p_FDR_corr_neg(:,:,k3) = beam.R.roi2voxel_tfm * beam.snpm.roi_p_duration_corr_neg(:,:,k3);
            end
        elseif isfield(beam.snpm,'p_uncorr_neg')
            [beam.snpm.p_FDR_corr_neg,beam.snpm.cutoff_neg] = nut_FDR(beam.snpm.p_uncorr_neg,q);
        end        
    end
    if isfield(beam,'ttest') && ~isfield(beam.ttest,'p_FDR_corr')        
        if isfield(beam.ttest,'FDR'), q=beam.ttest.FDR;
        else beam.ttest.FDR=q;
        end                 
        if isfield(beam.ttest,'roi_p_uncorr')
            [beam.ttest.roi_p_duration_corr,beam.ttest.cutoff] = nut_FDR(beam.ttest.roi_p_uncorr,q);
            beam.ttest.p_FDR_corr  = ones(size(beam.s{1}));
            for k3=1:size(beam.ttest.roi_p_uncorr,3)
                beam.ttest.p_FDR_corr(:,:,k3) = beam.R.roi2voxel_tfm * beam.ttest.roi_p_duration_corr(:,:,k3);
            end
        elseif isfield(beam.ttest,'p_uncorr')
            [beam.ttest.p_FDR_corr,beam.ttest.cutoff] = nut_FDR(beam.ttest.p_uncorr,q);
        end
    end    
    P = beam;
    return
end


% Make sorted vector
%-----------------------------------------------------------------------
dimp = size(uP);
[Ps,idx] = sort(uP(:));

% Set Benjamini & Yeuketeli cV for independence/PosRegDep case
%-----------------------------------------------------------------------
cV   = 1; 
S    = length(Ps);

% Calculate FDR inequality RHS
%-----------------------------------------------------------------------
Fi  = (1:S)'/S*q/cV;

% Find threshold
%-----------------------------------------------------------------------
I = max(find(Ps<=Fi));
if isempty(I)
  u = 0;
else
  u = Ps(I);
end

%-"Corrected" p-values
%-----------------------------------------------------------------------
Qs    = Ps*S./(1:S)'*cV;    

% -"Adjusted" p-values
%-----------------------------------------------------------------------
P     = zeros(size(Ps));

P(S) = Qs(S);
for i = S-1:-1:1
    P(i) = min(Qs(i),P(i+1));
end

% Bring back to original order and dimension
%-----------------------------------------------------------------------
P(idx)=P;
P = reshape(P,dimp);

if nargin>2         % legacy compatibility
    acc = (uP < u);
    oldout = P;
    P = acc;
end
