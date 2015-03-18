function beamout = fcm_basebeam(beam)

if ~isfield(beam,'baseline'), error('no baseline found.'), end

idx = (beam.timepts>=beam.baseline(1) & beam.timepts<beam.baseline(2));

if isfield(beam,'labels') && isfield(beam.labels,'contrasts'), beam.labels=rmfield(beam.labels,'contrasts'); end
if isfield(beam,'sinfo'), beam=rmfield(beam,'sinfo'); end

for s=1:length(beam.s)
    B = mean(beam.s{s}(:,idx,:),2);
    A = beam.s{s};
    
    beamout(s)=beam;
    beamout(s).s=[];
    beamout(s).s{1}=A;
    beamout(s).s{2}=repmat(B,[1 size(A,2) 1]);
    beamout(s).s{3}=zeros(size(A));
    
    if isfield(beam,'rois')
        B = mean(beam.rois{s}(:,idx,:),2);
        A = beam.rois{s};
        
        beamout(s).rois=[];
        beamout(s).rois{1}=A;
        beamout(s).rois{2}=repmat(B,[1 size(A,2) 1]);
        beamout(s).rois{3}=zeros(size(A));
    end

end



    