ic = 50;
clear temp data cfg

load testW_cog64
load data/chanlocs64.mat
temp = eeg_emptyset();
temp.icaweights = W;
temp.icasphere = sphere;
temp.icawinv = pinv(W*sphere);
temp.chanlocs = chanlocs;
temp.nbchan = length(W);
temp.icachansind = 1:temp.nbchan;
temp.dipfit.coord_transform = [];
data = eeglab2fieldtrip(temp,'componentanalysis','dipfit');
data.trial = {temp.icawinv(:,ic)};
data.time = {0};


load vis/cognionicsHeadModel.mat
load vis/cognionicssLORETA.mat
load(['vis/' cognionicsHeadModel.surfaces]);
vertices = surfData(3).vertices(setdiff(1:end,rmIndices),:);
nVert = size(vertices,1);
cfg.numdipoles = 1;
cfg.gridsearch = 'yes';
cfg.nonlinear  = 'yes';
cfg.grid.pos = vertices;
cfg.grid.inside = 1:nVert;
cfg.grid.leadfield = cell(nVert,1);
for it = 1:nVert
    cfg.grid.leadfield{it} = K(:,[it,it+nVert,it+2*nVert]);
end
cfg.model = 'moving';
cfg.latency = 'all';
cfg.component = ic;
cfg.channel = {chanlocs.labels};
cfg.vol.unit = 'mm';
cfg.vol.type = 'dipoli';
cfg.vol.skin_surface = 1;
cfg.vol.inner_skull_surface = 3;
cfg.vol.cond = [0.3300 0.0041 0.3300];
nVert_total = 0;
for it = 1:3
    cfg.vol.bnd(it).pnt = surfData(it).vertices;
    cfg.vol.bnd(it).tri = surfData(it).faces;
    nVert_total = nVert_total + length(surfData(it).vertices);
end
cfg.vol.mat = rand(size(W,1),nVert_total);