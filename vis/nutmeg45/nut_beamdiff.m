function nut_beamdiff(beamfile1,beamfile2,outfile)
% nut_beamdiff(beamfile1,beamfile2,outfile)
% saves output of beam1 - beam2 to outfile.
% note that beamfile1 and beamfile2 must contain the same voxels

ignoreorientation = true; % i can't imagine our sketchy orientations are useful for this

beam=load(beamfile1);
beam=nut_beam_legacy_compatibility(beam);
S1=beam.s; 
if isfield(beam,'rois'), R1=beam.rois; end

beam2=load(beamfile2);
beam2=nut_beam_legacy_compatibility(beam2);
S2=beam2.s; 
if isfield(beam2,'rois'), R2=beam2.rois; end

if ~isequal(beam.voxels,beam2.voxels)
    [beam.voxels,i1,i2] = intersect(beam.voxels,beam2.voxels,'rows');
    for k1=1:size(S1,1)
        for k2=1:size(S1,2)
            S1{k1,k2}=S1{k1,k2}(i1,:,:,:); 
            S2{k1,k2}=S2{k1,k2}(i2,:,:,:);
        end
    end    
end

if ndims(beam.s{1})>3 && ignoreorientation
    S1 = sqrt(sum(S1.^2,4));
    S2 = sqrt(sum(S2.^2,4));
%     s_beam1 = sqrt(s_th1.^2 + s_ph1.^2 + s_z1.^2);
%     s_beam2 = sqrt(s_th2.^2 + s_ph2.^2 + s_z2.^2);
end

if ~isequal(size(S1{1}),size(S2{1}))
    if isfield(beam,'rois'), error('unequal data size not supported for ROI data.'), end
    s1=size(S1{1}); s2=size(S2{1});
    if numel(s1)<numel(s2)
        s1(end+1:numel(s2))=1;
    elseif numel(s1)>numel(s2)
        s2(end+1:numel(s1))=1;
    end
    mult = s1./s2;
    if any(mult>1)
        if any(mult>1 & s2~=1), error('data dimensions do not match.'), end
        for k1=1:size(S2,1)
            for k2=1:size(S2,2)
                S2{k1,k2}=repmat(S2{k1,k2},mult);
            end
        end
    elseif any(mult<1)
        if any(mult<1 & s1~=1), error('data dimensions do not match.'), end
        for k1=1:size(S1,1)
            for k2=1:size(S1,2)
                S1{k1,k2}=repmat(S1{k1,k2},1./mult);
            end
        end        
    end
end

for k1=1:size(S1,1)
    for k2=1:size(S1,2)
        S1{k1,k2}=S1{k1,k2}-S2{k1,k2};
    end
end
beam.s=S1; clear S1 S2

if isfield(beam,'rois')
    for k1=1:size(R1,1)
        for k2=1:size(R1,2)
            R1{k1,k2}=R1{k1,k2}-R2{k1,k2};
        end
    end
    beam.rois=R1;
end

if isfield(beam,'labels')
    fn = fieldnames(beam.labels)';
    for k=find(~strcmpi(fn,'contrasts'))
        if ~any(strcmpi(strtok(beam.labels.(fn{k}),' '),{'Frequency' 'Time'}))
            beam.labels.(fn{k})=['\Delta ' beam.labels.(fn{k})];
        end
    end
end
    
if nargin<3
    [pathname,filename,fileext] = fileparts(beamfile1);
    outfile = fullfile(pathname,[filename '-diff.mat']);
end

save(outfile,'beam');