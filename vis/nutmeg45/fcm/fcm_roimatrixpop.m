function MM=fcm_roimatrixpop(voi,typ,doz)
% MM=fcm_roimatrixpop(voi,typ,do_z_tf)
%  typ      coherence type (has no effect for other functional connectivity measures)
%           1 = abs imag coh (default), 2 = mag squared coh, 3 = abs real coh
%  do_z_tf  transform connectivity values to z-scores (default is false).

if length(unique(voi.condnr))>1, error('Only one condition allowed'), end
if nargin<2 || isempty(typ), typ=1; end
if nargin<3, doz=false; end

cwd = pwd;
load(strrep(voi.pathnames{1},'s_beamtf_',''));
cd(fileparts(voi.pathnames{1}));
if isfield(CC,'roidef')
    load(CC.roidef);
else
    load(CC.roifile);
    load(R.roidef);
end

numroi = length(ROI.label);
numpat = length(voi.subjnr);
[dum,nt,nf] = size(CC.coh);
zo = double( ~( strcmp(CC.method,'ccohere') && typ==1 ) );
MM.coh = nan([numpat numroi numroi nt nf]);

for k=1:numpat
    if k~=1
        load(strrep(voi.pathnames{k},'s_beamtf_',''));
    end
    cd(fileparts(voi.pathnames{k}));
    load(CC.roifile)
    MM.coh(k,R.goodroi,R.goodroi,:,:) = fcm_roimatrix(CC,typ,doz,false);
end
MM.method = CC.method;
MM.frq = CC.frq;
MM.time = CC.time;
MM.roilabel = ROI.label;

cd(cwd);