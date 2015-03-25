function hmObj = makeHeadModel(name,chanlocs,output_dir)

if ~ischar(name)
    error('name must be a string'); end

% eeglabPath = fileparts(which('eeglab'));
% mobilabPath = [headModelPath filesep 'plugins' filesep 'mobilab'];
mobilabPath = fileparts(which('headModel'));

templateFile = [mobilabPath filesep 'data' filesep 'head_modelColin27_4825.mat'];
template = load([mobilabPath filesep 'data' filesep 'head_modelColin27_4825.mat']);
surfFile = [name '_Colin27_4825.mat'];
surfData = template.surfData;
save([output_dir filesep surfFile],'surfData');

%% display template head model: Colin27 + aal atlas + emotiv montage
dipfitdefs;
vol.bnd.pnt = surfData(3).vertices;
vol.bnd.tri = surfData(3).faces;
index = num2str(randi(1e5));
save(['temp' index],'vol')
rmpath(genpath([mobilabPath '/dependency/iso2mesh']))
electrodes = coregister(chanlocs,template_models(2).chanfile,'mesh',['temp' index]);
addpath(genpath([mobilabPath '/dependency/iso2mesh']))
delete(['temp' index '.mat'])

% elec = 1000*[-[chanlocs.Y]; [chanlocs.X]; [chanlocs.Z]]';
elec = electrodes.pnt;
label = {chanlocs.labels};
hmObj = headModel('surfaces',surfFile,'atlas',template.atlas,'fiducials',template.fiducials,'channelSpace',elec,'label',label);
% plotHeadModel(hmObj);

%% warping channel space to template
aff = hmObj.warpChannelSpace2Template(templateFile,surfFile,'affine');
% plotHeadModel(hmObj);

%% solving the forward problem with OpenMEEG

dir_original = cd;
cd(output_dir);

conductivity = [0.33 0.022 0.33]; % brain and scalp = 0.33 S/m, skull = 0.022 S/m; these conductivies were taken from
                                  % Valdes-Hernandez et al., 2009, Oostendrop TF, 2000; Wendel and Malmivuo, 2006      
normal2surface = true;
hmObj.computeLeadFieldBEM(conductivity,normal2surface);

cd(dir_original);

%% solving the inverse problem with sLORETA
% K: lead field matrix
% L: Laplaciian operator
% rmIndices: indices to be removed (the Thalamus)
% surfData(3).vertices: source space
addpath(genpath('/home/lpiontonachini/Desktop/eeglab/plugins/mobilab/dependency/'))

[sourceSpace,K,L,rmIndices] = getSourceSpace4PEB(hmObj);

%--
[U,S,V] = svd(K/L,'econ');
Ut = U';
s2 = diag(S).^2;
iLV = L\V;

eval([name '_HeadModel = hmObj;']);
save([output_dir filesep name '_HeadModel'],[name 'HeadModel'])
save([output_dir filesep name '_LFMetc'],'K','L','rmIndices','Ut','s2','iLV')




