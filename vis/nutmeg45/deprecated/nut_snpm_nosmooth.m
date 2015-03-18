function [Tc,MaxTc,nPc]=nut_snpm_nosmooth(X,Q,VOX,SP)
% NUT_SNPM_NOSMOOTH  performs SnPM permutations
%
%    [T,MaxT,nP]=nut_snpm(X,Q,VOX,SP)
%
% X         spatially normalized data (subject*voxels)
% Q         vector of good voxels (present in all subjects and non-0 mean)
% VOX       voxelsize in mm [x y z]
% SP        SnPM parameter structure containing the following fields:
%           CONT      Condition parameter
%           STAT      'T' if t-tests, 'F' if ANOVA
%           PiCond    Permutation matrix
%           sHCform   function to update design matrix with spm_DesMtx
%           numperm   number of actual permutions
%           hPerms    true if permutation matrix is symmetric.
%           iSubj     subject labeling vector (ANOVA only)
%           iRepl     condition labeling vector (ANOVA only)
%
% This function is called by nut_timef_stats (GUI) and nut_beamstats.
% Adapted from the function snpm_cp in the SnPM3 toolbox:
%    http://www.fil.ion.ucl.ac.uk/spm/snpm

error('This function has been integrated into nut_snpm and is obsolete.')

% Initialize
%-----------
perm = 1;
HC = eval(SP.sHCform);
df   = size(HC,1) - rank(HC);
numgvox = size(X,2);

nPc   = ones(numgvox,1); 
%Tc    = zeros(numgvox,1);
MaxTc = zeros(SP.numperm,2);       % permutations*[max min]

% Variance
%---------
fprintf('\tCalculating variance of each voxel...\n')
BETA  = pinv(HC)*X;      % X: numsubj*num_non0voxels
ResSS = sum((X - HC*BETA).^2);

% Calculate t- or F-values
%--------------------------
fprintf('\tCalculating %s-values...\n',SP.STAT)
Co = SP.CONT;

if SP.STAT=='T'
    Tc(:,1) = Co*BETA ./ sqrt((ResSS*(Co*pinv(HC'*HC)*Co'))/df);
else        % F!
    %pX   = pinv(HC);
	Tc(:,1) = (sum(((Co*BETA)'*inv(Co*pinv(HC'*HC)*Co'))' .* ...
			(Co*BETA),1)/rank(Co)) ./ (ResSS/df);
end

MaxTc(1,:) = [max(Tc) -min(Tc)];


% Permutations
%-------------  
fprintf('\tPerforming permutations...\n')
%TP   = zeros(1,size(BETA,2));
for perm = 2:SP.numperm
    % Calculate variance
    HC = eval(SP.sHCform);
    BETA  = pinv(HC)*X;
    ResSS = sum((X - HC*BETA).^2);

    % Calculate t- or F-values
    if SP.STAT=='T'
        TPc(:,1) = Co*BETA ./ sqrt((ResSS*(Co*pinv(HC'*HC)*Co'))/df);
    else        % F!
        %pX   = pinv(HC);
	    TPc(:,1) = (sum(((Co*BETA)'*inv(Co*pinv(HC'*HC)*Co'))' .* ...
			    (Co*BETA),1)/size(Co,1)) ./ (ResSS/df);
    end
    MaxTc(perm,:) = [ max(TPc) -min(TPc) ];
    %MaxT(perm,1) = max(TP(1,:));
    
    if SP.hPerms
        nPc = nPc + (TPc>=Tc) + (-TPc>=Tc);   % NB: Worry if T=TP=0
        % if STAT=='T', then TP, T >=0, so (-TP>=T) will be empty.
    else
        nPc = nPc + (TPc>=Tc);
    end

end             


 