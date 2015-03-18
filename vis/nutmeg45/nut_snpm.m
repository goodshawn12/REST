function [T,p_uncorr,MaxT,MaxCluster,th1]=nut_snpm(X,Q,voxels,voxelsize,SP)
% NUT_SNPM  performs SnPM permutations
%
%    [perm1,permMax]=nut_snpm(X,Q,voxels,voxelsize,SP)
%
% X         spatially normalized data (subject*voxels)
% Q         vector of good voxels (present in all subjects and non-0 mean)
% voxels    MNI coordinates in mm of common voxels
% voxelsize voxelsize in mm [x y z]
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


% Initialize
%-----------
perm = 1;
HC = eval(SP.sHCform);
G  = SP.G; 
df = size(HC,1) - rank([HC G]);

% Ensure validity of contrast of condition effects, zero pad
if size(SP.CONT,1)==1
    if size(HC,2)>1
        SP.CONT(1:size(HC,2)) = SP.CONT(1:size(HC,2)) - mean(SP.CONT(1:size(HC,2)));
    end
end
% Zero pad for G partitions.
SP.CONT  = [SP.CONT zeros(1,size(G,2))];
% Kill zero-rank components
if rank(SP.CONT)<size(SP.CONT,1)
  [u s] = spm_svd(SP.CONT'); 
  SP.CONT = full(u*sqrt(s))'; clear u s
end

numgvox = size(X,2);
if (SP.dosmoothvar && length(Q)~=numgvox), error('Wrong Q'), end
voxelsblob = mm2voxels(voxels,voxelsize);
dimxyz = max(voxelsblob);
voxelsblob = voxelsblob(Q,:);

docluster = ( max(max(abs(round(voxelsblob)-voxelsblob))) < 1e-2 );
if ~docluster, warning('Non-integer voxel coordinates. Cluster correction cannot be calculated.')
else docluster = (numgvox==length(Q));  % No cluster correction for ROIs.
end

p_uncorr   = ones(numgvox,1);
MaxT = zeros(SP.numperm,(SP.STAT~='F')+1);       % permutations*[max min]
MaxCluster = zeros(SP.numperm,(SP.STAT~='F')+1);

% Main Part
%----------
fprintf('\tWorking across subjects...\n')
if SP.STAT=='R'
    goodvox = ( sum( X==0 , 1) < 3 );
    for vv=find(goodvox)  %1:numgvox
        X(:,vv) = X(:,vv) - mean(X(:,vv));
        X(:,vv) = X(:,vv) ./ norm(X(:,vv));
    end
    X(:,~goodvox)=0;
end
BETA  = pinv([HC G])*X;      % X: numsubj*num_non0voxels

if SP.STAT~='R'
    ResSS = sum((X - [HC G]*BETA).^2);
end
if SP.dosmoothvar
    fprintf('\tSmoothing variance with a 20x20x20 mm gaussian kernel...\n')     

    vFWHM     = [20 20 20];
    SmResSS   = zeros(dimxyz);
    SmMask    = zeros(dimxyz);
    TmpVol    = zeros(dimxyz);
    TmpVol(Q) = ones(size(Q));  % Q: non0voxels
    spm_smooth(TmpVol,SmMask,vFWHM./voxelsize);
    TmpVol(Q) = ResSS;
    spm_smooth(TmpVol,SmResSS,vFWHM./voxelsize);
    ResSS     = SmResSS(Q)./SmMask(Q);
end

% Calculate t- R- or F-values
%--------------------------
fprintf('\tCalculating %s-values...\n',SP.STAT)
Co = SP.CONT;

switch SP.STAT
    case 'T'
        T(:,1) = Co*BETA ./ sqrt((ResSS*(Co*pinv([HC G]'*[HC G])*Co'))/df);
        MaxT(1,:) = [max(T) -min(T)];    
    case 'F'
        %pX   = pinv(HC);
        T(:,1) = (sum(((Co*BETA)'*inv(Co*pinv([HC G]'*[HC G])*Co'))' .* ...
                (Co*BETA),1)/rank(Co)) ./ (ResSS/df);
        MaxT(1) = max(T);
    case 'R'
        T(:,1) = BETA;
        MaxT(1,:) = [max(T) -min(T)];
    otherwise
        error('I do not know the stats you requested!?...')
end

if SP.hPerms
    p_uncorr = p_uncorr + (T<0);
end

if docluster
    Ts = sort(T);
    th1= Ts(round(0.95*length(Ts)));

    firstpass = (T > th1);
    if any(firstpass)
        MaxCluster(1,1) = max(spm_max(T(firstpass).',voxelsblob(firstpass,:).'));
    end
    if SP.STAT~='F'
        firstpass = (T < -th1);
        if any(firstpass)
            MaxCluster(1,2) = max(spm_max(-T(firstpass).',voxelsblob(firstpass,:).'));
        end
    end
else
    MaxCluster=nan(SP.numperm,size(MaxT,2));
    th1=[];
end

% Permutations
%-------------  
fprintf('\tPerforming permutations...\n')
for perm = 2:SP.numperm
    % Calculate variance
    HC = eval(SP.sHCform);
    BETA  = pinv([HC G])*X;
    if SP.STAT~='R'
        ResSS = sum((X - [HC G]*BETA).^2);
    end
    if SP.dosmoothvar
        % Smooth Variance
        TmpVol(Q) = ResSS;
        spm_smooth(TmpVol,SmResSS,vFWHM./voxelsize);
        ResSS     = SmResSS(Q)./SmMask(Q);
    end

    % Calculate t- R- or F-values
    if SP.STAT=='T'
        TP(:,1) = Co*BETA ./ sqrt((ResSS*(Co*pinv([HC G]'*[HC G])*Co'))/df);
        MaxT(perm,:) = [ max(TP) -min(TP) ];
    elseif SP.STAT=='F' 
        %pX   = pinv(HC);
	    TP(:,1) = (sum(((Co*BETA)'*inv(Co*pinv([HC G]'*[HC G])*Co'))' .* ...
			    (Co*BETA),1)/size(Co,1)) ./ (ResSS/df);
        MaxT(perm) = max(TP);
    else
        TP(:,1) = BETA;
        MaxT(perm,:) = [ max(TP) -min(TP) ];
    end
    
    if SP.hPerms
        p_uncorr = p_uncorr + (TP>=T) + (-TP>=T);   % NB: Worry if T=TP=0
        % if STAT=='T', then TP, T >=0, so (-TP>=T) will be empty.
    else
        p_uncorr = p_uncorr + (TP>=T);
    end
    
    if docluster
        firstpass = (TP > th1);
        if any(firstpass)
            MaxCluster(perm,1) = max(spm_max(TP(firstpass).',voxelsblob(firstpass,:).'));
        end
        if SP.STAT~='F'
            firstpass = (TP < -th1);
            if any(firstpass)
                MaxCluster(perm,2) = max(spm_max(-TP(firstpass).',voxelsblob(firstpass,:).'));
            end
        end    
    end
end             

%===========
function voxelsblob = mm2voxels(voxelsMNI,voxelsize)

voxelsMNI=double(voxelsMNI);

 blob2mri_tfm = [ voxelsize(1)        0                 0             min(voxelsMNI(:,1))-voxelsize(1)
                            0          voxelsize(2)        0             min(voxelsMNI(:,2))-voxelsize(2)
                            0                 0         voxelsize(3)     min(voxelsMNI(:,3))-voxelsize(3)
                            0                 0                 0                                1 ];
voxelsblob = nut_coordtfm(voxelsMNI,inv(blob2mri_tfm));
 