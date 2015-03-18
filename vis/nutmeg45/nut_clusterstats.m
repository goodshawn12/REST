function CP = nut_clusterstats(R,P,beam,MinClusterSize,firstpass)
% nut_clusterstats  thresholds multiple voxel stats with cluster size.
%
% p_cluster_corr = nut_clusterstats(statval,p_uncorr,beam,MinClusterSize,fistpass)
%    statval          e.g., T values or correlation coefficients R
%    p_uncorr         p-values not corrected for multiple testing
%    beam             NUTMEG beam (or nuts) structure with voxel coordinates and voxelsize
%    MinClusterSize   Minimum cluster size of significant voxels (scalar or 
%                     numtimepts*numfreqbands matrix).
%    firstpass        (optional) threshold for cluster formation. Scalar with min/max statval or
%                     logical array of same size as statval. If no value is
%                     indicated, firstpass is (p_uncorr < 0.05)
% 
% Calls spm_max and spm_clusters

sigPu = (P<0.05);
sigRu = R(sigPu);
ispos = any(sigRu(:)>0);
isneg = any(sigRu(:)<0);
isc   = isscalar(MinClusterSize);

if nargin<5 || isempty(firstpass)
    firstpass = sigPu;
elseif isscalar(firstpass)
    firstpass = (abs(R)>abs(firstpass));
end

[numvox,numtim,numfrq] = size(R);
voxelsblob = nut_voxels2blob(beam);
if max(max(round(voxelsblob(1:5,:))-voxelsblob(1:5,:)))>.05, error('Cannot do clusterstats with non-integer coordinates.'), end
voxelsblob=round(voxelsblob);
CP = P;
CP(~firstpass) = 1;

for ff = 1:numfrq
    for tt=1:numtim
        firstpassC = find(firstpass(:,tt,ff));
        if ~isempty(firstpassC)
            X = R(firstpassC,tt,ff).';
            XYZ = voxelsblob(firstpassC,:).';
            clusterpass = false(size(X));
            if ispos
                ixpos=(X>0);
                if any(ixpos)
                    [N dum dum AR]=spm_max(X(ixpos),XYZ(:,ixpos));
                    AV = spm_clusters(XYZ(:,ixpos));
                    clusterpass(ixpos) = ismember(AV,unique(AR( N >= MinClusterSize(tt-isc*(tt-1),ff-isc*(ff-1)) )));
                end
            end
            if isneg
                ixneg=(X<0);
                if any(ixneg)
                    [N dum dum AR]=spm_max(-X(ixneg),XYZ(:,ixneg));
                    AV = spm_clusters(XYZ(:,ixneg));
                    clusterpass(ixneg) = ismember(AV,unique(AR( N >= MinClusterSize(tt-isc*(tt-1),ff-isc*(ff-1)) )));
                end
            end
            CP(firstpassC(~clusterpass),tt,ff) = 1;
        end
    end
end